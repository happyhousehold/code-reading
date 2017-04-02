% this file implements an APL method for solving  %
% 1/2 * ||A*x - b||^2 

clear;
%  load n4000m3000r1
% load n8000m4000r1;
% load n8000m4000gaussian;
% load n4000m3000gaussian;
% % load Bdata2062x2062_1;
% load Bdata2062x4124_3;
% data.A=Bdata.A;
% data.b=Bdata.rhs;
% test_x=Bdata.sol;
    
data.A =randn(3000,10000);
data.R =1.1;
test_x=randn(10000,1);
 test_x=test_x/norm(test_x);
data.b= data.A * test_x;
data.n=size(test_x,1);

% data.b= data.A * test_x + 1 * randn(8000,1);
% data.b= randn(4000,1)+1e-4 * randn(4000);
% data.A=n4000m3000gaussian.data.A;
% data.b=n4000m3000gaussian.data.b;
% test_x=n4000m3000gaussian.data.sol;
% data.n = size(n4000m3000gaussian.data.sol, 1);


% setting the control parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% change theta, beta to 0.2 & 0.1 to update LB
control.theta = 0.2;
control.beta = 0.2;
control.epsilon = 1e-15;

control.iter_limit = 3000;
max_inner_iter=control.iter_limit;
control.bundle_limit =5;
% control.M1 = perm_gen (control.bundle_limit+1); 
control.M2 = perm_gen (control.bundle_limit+2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% initialization
% compute the initial lower and upper bound

p0 = zeros(data.n, 1);

% p0 = randn(data.n, 1);
% p0=p0/norm(p0)*data.R;

[f,g] = oracleqp(data, p0);
const = f - g' * p0;
% compute t_p0 \in argmin{ f(p_0) + g(p_0)' * (x - p0): \|x\| \le R }
t_p0  = - data.R * g / norm(g);


[ft,gt] =  oracleqp(data, t_p0 );
 LB   = f + g' * (t_p0 - p0);
% LB=0;
UB   = min(f,ft);



% select the initial prox center
if (ft < f),
    x_ubt = t_p0;
else
    x_ubt = p0;
end;
prox_center = x_ubt;
nstep = 0;
s = 1;

start_APL=clock;
Bundle.size = 0;
Bundle.matrix = zeros(control.bundle_limit, data.n);
Bundle.const = zeros(control.bundle_limit, 1);

terminate = 0;

%%%%%%%%%%%%%%%%%%%  step 0  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (terminate==0)

    if (UB - LB <= control.epsilon || nstep > control.iter_limit )
        break;
    end
    x_ubt = prox_center;
    x_t = prox_center;
    LS = control.beta * LB + (1 - control.beta) * UB;
    stepLB = LB;
    stepUB = UB;
    
    % set up the constrain < x_t - prox_center, x_t - x >  <= 0
    BarX.a = zeros(data.n, 1);
    BarX.b = 0;
    
    t=0;
    
    %%%%%%%%%%%%%%%%%%%%%   step 1   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inner_iter=0;
    while 1 == 1
        inner_iter=inner_iter+1;
        nstep = nstep + 1;
        t     = t + 1;
        alpha_t = 2.0 / (t + 1.0);
        
%         if (stepUB <= control.epsilon || nstep > control.iter_limit )
        if (stepUB-stepLB <= control.epsilon || nstep > control.iter_limit )
            terminate = 1;
            break;
        end
        %%% updata the lower bound
        
        x_lbt = (1 - alpha_t) * x_ubt + alpha_t * x_t;
        [f,g]= oracleqp(data, x_lbt);
        
        newg = g';
        newconst = f - g' * x_lbt;
        
        fprintf('nstep=%d, UB=%.2e, LB=%.2e, stepUB=%.2e, stepLB=%.2e, GAP=%.2e, Bundle=%d\n',...
            nstep, UB, LB, stepUB, stepLB, stepUB - stepLB, Bundle.size);

        
        %%%%%%%%%%%%%%%%%%%%  step  2    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% compute prox-mapping
        
        % using CVX
        %      x_t = ProxMappingCVX(domain, Bundle, data, BarX, c, LS);
        % using IPM
        %      x_t = ProxMapping(data, Bundle, control, BarX, prox_center, LS);
        
        %%%% using kkt_solver
         [x_t,er] = ProxMapping_kkt (x_t, Bundle, BarX, newg, newconst, prox_center, LS, control);
        
       if  (er==1 || norm(x_t)> data.R) % data.R should change to  2^0.5 * data.R
            LB = LS ;     
            UB = stepUB;
            prox_center= x_ubt;
            % reset Bundle here is necessary
            Bundle.size = 0;
            Bundle.matrix = zeros(control.bundle_limit, data.n);
            Bundle.const = zeros(control.bundle_limit, 1);
            break;
       end

        %%%%%%%%%%%%%%%%%%%%  step  3    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmp_xub = alpha_t * x_t + (1 - alpha_t) * x_ubt;
        [fu, gu] = oracleqp(data, tmp_xub);
        
        % update the upper bound
        if (fu < stepUB) 
            stepUB = fu;
            x_ubt  = tmp_xub;

        end
        
        
        % update the Bundle
         Bundle.matrix = [newg ; Bundle.matrix(1 : control.bundle_limit-1, :)];
         Bundle.const = [newconst ; Bundle.const(1 : control.bundle_limit-1)];
         if Bundle.size <  control.bundle_limit
            Bundle.size = Bundle.size + 1;
         end
        
        % update the constraint < BarX.a, x > <= BarX.b
        BarX.a = prox_center - x_t;
        BarX.b = BarX.a' * x_t;
        
        % significant progress on the upper bound
        if stepUB <= LS + control.theta * (UB - LS)
            LB = stepLB;
            UB = stepUB;
            prox_center  = x_ubt;
            break;
        end
        
        if inner_iter==max_inner_iter
            LB = stepLB;
            UB = stepUB;
            prox_center  = x_ubt;
            break;
        end
              
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    end
    
end


 
time_APL=etime(clock,start_APL);
fprintf('end of execution.\n');
str  = sprintf('nstep=%d, UB=%.2e, LB=%.2e, stepUB=%.2e, stepLB=%.2e, time=%5.2f\n',...
    nstep, UB, LB, stepUB, stepLB, time_APL);
disp(str);



