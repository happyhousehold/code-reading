%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MATLAB Code for the Accelerated Prox-Level (APL) algorithm            %
%     Author: Guanghui (George) Lan                                      %
%     Institute: University of Florida, Industrial & Systems Engineering %
%     @All rights reserved 2010                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this file implements an APL algorithm for solving QP %
function APL_Eigen_basic(control, domain, data)

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
s = 0;

start_APL=clock;
% define the bundle structure
Bundle.size = 0;
Bundle.matrix = zeros(control.bundle_limit, domain.n);
Bundle.const = zeros(control.bundle_limit, 1);
lb_dual = 0;
map_dual = 0;

basemu = control.theta * (1-control.lambda) / (2*log(data.m));

terminate = 0;
while terminate == 0 
    if (UB - LB <= control.epsilon || nstep > control.iter_limit )
        break;
    end
    
    % start stage s = 1, 2,..
    s = s+1;
    xt = x_ubt;
    c  = x_ubt;
    LS = LB + control.lambda * (UB - LB);
    stepLB = LB;
    stepUB = UB;
    
    % update the smoothing parameter
    data.mu = basemu * (UB-LB);

%    if data.mu < control.epsilon
%        data.mu = control.epsilon;
%    end
    
 %   data.mu = 0;
    % clear up the constraint that defines \bar(X)
    BarX.a = zeros(domain.n, 1);
    BarX.b = 0;
       
    t = 0;

    while 1 == 1   
        nstep = nstep + 1;
        t     = t + 1;
        alpha_t = 2.0 / (t + 1.0);
         
        if (stepUB - stepLB <= control.epsilon || nstep > control.iter_limit )
            terminate = 1;
            break;
        end
        
        % obtain an initial feasible solution to X_t for the computation of
        % the lower bound
       
        x_lbt = (1 - alpha_t) * x_ubt + alpha_t * xt;
        [f,g,f0,dlb]= oracle(data, x_lbt, 0);
        
        disp(sprintf('nstep=%d, UB=%.6e, LB=%.6e, stepUB=%.6e, stepLB=%.6e, Gap=%.6e, Bundle=%d\n',...
          nstep, UB, LB, stepUB, stepLB, stepUB-stepLB,Bundle.size));
        if nstep > control.iter_limit
            break;
        end

        if data.DomainType == 0,
            hstar = ComputeNewLB(domain, Bundle, data, BarX, f, g, x_lbt, LS, control.lb_mode);
        else,
            hstar = ComputeNewLB_box(domain, Bundle, data, BarX, f, g, x_lbt, LS, control.lb_mode);
        end;
  %      hstar = ComputeNewLBCVX(domain, Bundle, data, BarX, f, g, x_lbt, LS, control.lb_mode);
        
        stepLB = max(stepLB, min(hstar,LS));
 
  %      if stepLB <dlb,
  %          stepLB = dlb;
  %      end

        
  %      if control.bundle_limit > 0,                     
  %              if Bundle.size < control.bundle_limit
  %                      Bundle.size = Bundle.size + 1;
  %                      Bundle.matrix(Bundle.size, :) = g';
  %                      Bundle.const(Bundle.size) = f - g'* x_lbt;
  %                  else
  %                      Bundle.matrix(1:Bundle.size-1, :) = Bundle.matrix(2:Bundle.size, :);
  %                      Bundle.const(1:Bundle.size-1) = Bundle.const(2:Bundle.size);
  %                      Bundle.matrix(Bundle.size, :) = g';
  %                      Bundle.const(Bundle.size) = f - g'* x_lbt;
  %                  end
  %      end
        
        if control.bundle_limit > 0,                     
                if Bundle.size < control.bundle_limit
                        Bundle.size = Bundle.size + 1;
                        newguy = Bundle.size;                        
                else  % remove the oldest guy
                        newguy = mod(nstep-1, control.bundle_limit)+1;
                end
                Bundle.matrix(newguy, :) = g';
                Bundle.const(newguy) = f - g'* x_lbt;
        end
        
        % significant progress on the lower bound
        if stepLB >= LS - control.theta * (LS - LB),
            LB = stepLB;
            UB = stepUB;
            c  = x_ubt;
            break;
        else,
        %    xt = ProxMappingCVX(domain, Bundle, data, BarX, f, g, x_lbt, c, LS);
            if data.DomainType == 0,
                xt = ComputeProxMapping(domain, Bundle, data, BarX, f, g, x_lbt, c, LS);
            else
                if data.DomainType == 1,
                    xt = ComputeProxMapping_box(domain, Bundle, data, BarX, f, g, x_lbt, c, LS);
                else
                    xt = ComputeProxMapping_box_entropy(domain, Bundle, data, BarX, f, g, x_lbt, c, LS);
                end
            end;
  
            tmp_xub = alpha_t * xt + (1 - alpha_t) * x_ubt;
            [fu, gu, f0u] = oracle(data, tmp_xub, 1);
            if f0u < stepUB,
                stepUB = f0u;
                x_ubt = tmp_xub;
            end
            
            % the constraint that defines \bar{X}
            % BarX.a = xt - c;
            % the constraint that defines \bar{X}
            % domega_c x_t (x - x_t) \ge 0
            % note that domega x_t = 
            if data.DomainType == 0,
                BarX.a = log(xt) - log(c);          
            else 
                if data.DomainType == 1,
                    BarX.a = xt - c;
                else
                    BarX.a = (log(xt/domain.n) - log(c/domain.n)) / domain.n;
                end
            end;
            
            BarX.b = BarX.a' * xt;
        
            
            % significant progress on the upper bound
            if stepUB <= LS + control.theta * (UB - LS),
                LB = stepLB;
                UB = stepUB;
                c  = x_ubt;
                break;
            end
        end            
    end
end
time_APL=etime(clock,start_APL);
disp(sprintf('end of execution.\n'));
str  = sprintf('nstep=%d, UB=%.6e, LB=%.6e, stepUB=%.6e, stepLB=%.6e, time=%5.2f\n',...
    nstep, UB, LB, stepUB, stepLB, time_APL);
disp(str);
frep = fopen('report.txt','a');
fprintf(frep,'Summary: %s**\n',str);
fclose(frep);




