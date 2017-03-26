%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MATLAB Code for the Accelerated Prox-Level (APL) algorithm            %
%     Author: Guanghui (George) Lan                                      %
%     Institute: University of Florida, Industrial & Systems Engineering %
%     @All rights reserved 2010                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this file implements an NERML algorithm for solving eigenvalue problem %
function NERML_Eigen(control, domain, data)

start_NERML=clock;
% compute the initial lower and upper bound
data.mu = 0;
if data.DomainType == 0, %simplex
    p0 = ones(domain.n,1)/domain.n;
    [f,g,f0] = oracle(data, p0, 0);
    UB   = f0; % initial upper bound
    % compute the initial lower bound
    % min f0(p_0) + <f0'(p_0), x - p_0>: x in X
    % X is a standard simplex < Bnd
    Bnd = 1;

    [v,I] = min(g);
    t_p0 = zeros(domain.n,1);
    t_p0(I) = Bnd;

    LB   = f0 + g' * (t_p0 - p0); % initial lower bound
else
    if data.DomainType == 1 %box [0,1] domain, and using \|.\|^2/2 prox-function
        p0 = zeros(domain.n,1);
        [f,g,f0] = oracle(data, p0, 0);
        UB = f0;
        % min f0(p_0) + <f0'(p_0), x - p_0>: x in X, X is box
        LB = f0 - g'*p0 + sum(min(0,g));
    else %box [0,1] domain, but we embed it into a 1 \le \sum x_i \le n using 
        % the prox-function \sum_i x/n log x/n
        p0 = ones(domain.n,1);
        [f,g,f0] = oracle(data, p0, 0);
        UB   = f0;
        % min f0(p_0) + <f0'(p_0), x - p_0>: x in X, X is box
        LB = f0 - g'*p0 + sum(min(0,g));
    end
end
% select the initial prox center
prox_center = p0;
x_ubt = p0;
nstep = 0;
s = 1;

% define the bundle structure
Bundle.size = 0;
Bundle.matrix = zeros(control.bundle_limit, domain.n);
Bundle.const = zeros(control.bundle_limit, 1);

while 1==1
    if (UB - LB <= control.epsilon || nstep > control.iter_limit )
        break;
    end
    
    % start stage s = 1, 2,..
    t = 1;
    xt = x_ubt;
    c = x_ubt;
    LS = LB + control.lambda * (UB - LB);
    stepLB = LB;
    stepUB = UB;
    
    % clear up the constraint that defines \bar(X)
    BarX.a = zeros(domain.n, 1);
    BarX.b = 0;
    t = 0;

    while 1 == 1   
        nstep = nstep + 1;


        [f,g,f0]= oracle(data, xt, 0); 
        if (f < stepUB),
            stepUB = f;
            x_ubt  = xt;
        end;
        
         % update bundle
        if control.bundle_limit > 0,                     
             if Bundle.size < control.bundle_limit
                  Bundle.size = Bundle.size + 1;
                  Bundle.matrix(Bundle.size, :) = g';
                  Bundle.const(Bundle.size) = f - g'* xt;
             else
                  Bundle.matrix(1:Bundle.size-1, :) = Bundle.matrix(2:Bundle.size, :);
                  Bundle.const(1:Bundle.size-1) = Bundle.const(2:Bundle.size);
                  Bundle.matrix(Bundle.size, :) = g';
                  Bundle.const(Bundle.size) = f - g'* xt;
            end
        end
            
        disp(sprintf('nstep=%d, UB=%.6e, LB=%.6e, stepUB=%.6e, stepLB=%.6e, gap=%.6e, Bundle=%d\n', ...
            nstep, UB, LB, stepUB, stepLB, stepUB - stepLB, Bundle.size));
        if nstep > control.iter_limit
            break;
        end
        % compute h* = min{ f + g' * (x - x_lb): x \in X_{t-1}      
        %hstar = ComputeNewLB(domain, Bundle, BarX, f, g, xt, LS, control.lb_mode);
        %hstar = ComputeNewLBCVX(domain, Bundle, data, BarX, f, g, xt, LS, control.lb_mode);
        if data.DomainType == 0,
            hstar = ComputeNewLB(domain, Bundle, data, BarX, f, g, xt, LS, control.lb_mode);
        else,
            hstar = ComputeNewLB_box(domain, Bundle, data, BarX, f, g, xt, LS, control.lb_mode);
        end;
        
        stepLB = max(stepLB, min(hstar,LS));
        
        % significant progress on the lower bound
        if stepLB >= LS - control.theta * (LS - LB),
            LB = stepLB;
            UB = stepUB;
            c  = x_ubt;
            break;
        else,            
            % compute prox-mapping
           % xtplus = ProxMappingCVX(domain, Bundle, data, BarX, f, g, xt,c, LS);
           if data.DomainType == 0,
                xtplus = ComputeProxMapping(domain, Bundle, data, BarX, f, g, xt, c, LS);
            else
                xtplus = ComputeProxMapping_box(domain, Bundle, data, BarX, f, g, xt, c, LS);
            end;
            % significant progress on the upper bound
            if stepUB <= LS + control.theta * (UB - LS),
                LB = stepLB;
                UB = stepUB;
                c  = x_ubt;
                break;
            end

           
            % the constraint that defines \bar{X}
            if data.DomainType == 0
                BarX.a = log(xt) - log(c);
            else
                BarX.a = xt - c;
            end
            BarX.b = BarX.a' * xt;
            
            xt = xtplus;
        end       
    end
end
time_NERML=etime(clock,start_NERML);
disp(sprintf('end of execution.\n'));
str  = sprintf('nstep=%d, UB=%.6f, LB=%.6f, stepUB=%.6f, stepLB=%.6f, time=%5.2f\n',...
    nstep, UB, LB, stepUB, stepLB, time_NERML);
disp(str);
frep = fopen('report.txt','a');
fprintf(frep,'Summary: %s**\n',str);
fclose(frep);




