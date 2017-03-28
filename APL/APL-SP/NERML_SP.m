%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MATLAB Code for the Accelerated Prox-Level (APL) algorithm            %
%     Author: Guanghui (George) Lan                                      %
%     Institute: University of Florida, Industrial & Systems Engineering %
%     @All rights reserved 2010                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this file implements an NERML algorithm for solving eigenvalue problem %
function NERML_SP(control, domain, data)

start_NERML=clock;
% define the bundle structure
Bundle.size = 0;
Bundle.matrix = zeros(control.bundle_limit, domain.n);
Bundle.const = zeros(control.bundle_limit, 1);
% clear up the constraint that defines \bar(X)
BarX.a = zeros(domain.n, 1);
BarX.b = 0;

% compute the initial lower and upper bound
data.mu = 0;
p0 = zeros(domain.n,1);
%%% compute the projection to the feasible set
p0 = ComputeProxMapping(domain, Bundle, data, BarX, 0, 0, 0, p0, -inf);
%p0 = data.iniSol;
[f,g,f0] = oracle(data, p0);

% compute the initial lower bound
LB = ComputeNewLB(domain, Bundle, data, BarX, f, g, p0, -inf, control.lb_mode);
UB = f0;

% select the initial prox center
prox_center = p0;
x_ubt = p0;
nstep = 0;
s = 1;



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


        [f,g,f0]= oracle(data, xt); 
        
        if (f < stepLB)
            disp(sprintf('UB<LB, possible numerical error!'));            
        else
            if (f < stepUB),
                stepUB = f;
                x_ubt  = xt;
            end
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
        hstar = ComputeNewLB(domain, Bundle, data, BarX, f, g, xt, LS, control.lb_mode);
        
        if min(hstar, LS) > UB,
            disp(sprintf('UB<LB, possible numerical error!'));
        else    
            stepLB = max(stepLB, min(hstar,LS));
        end
        
        % significant progress on the lower bound
        if stepLB >= LS - control.theta * (LS - LB),
            LB = stepLB;
            UB = stepUB;
            c  = x_ubt;
            break;
        else,            
            % compute prox-mapping
            
           xtplus = ComputeProxMapping(domain, Bundle, data, BarX, f, g, xt, c, LS);

            % significant progress on the upper bound
            if stepUB <= LS + control.theta * (UB - LS),
                LB = stepLB;
                UB = stepUB;
                c  = x_ubt;
                break;
            end
          
            % the constraint that defines \bar{X}
            BarX.a = xt - c;
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




