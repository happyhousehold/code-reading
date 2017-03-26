function [lb_t] = ComputeNewLB_box(domain, Bundle, data, BarX, f, g, x_lb, LS, mode)
clear lb_prob;
clear lb_res
clear lb_sol;

if Bundle.size == 0,
%        cvx_begin
%        cvx_quiet(true);
%            variable x(domain.n);
%            minimize (f + g' * (x - x_lb));
%            subject to
%                x == box(domain.n);
%                BarX.a' * x >= BarX.b;
%        cvx_end
     %specify the linear constrain matrix
    lb_prob.c = g;
    lb_prob.a = [-sparse(BarX.a')];
    lb_prob.blc = [-inf];
    lb_prob.buc = [-BarX.b];
    lb_prob.blx = zeros(domain.n,1); 
    %lb_prob.blx = -data.DomainType*ones(domain.n,1);
    lb_prob.bux = data.DomainType*ones(domain.n,1);
else,
    if mode == 0,
%            cvx_begin
%            cvx_quiet(true);
%                variable x(domain.n);
%                minimize (f + g' * (x - x_lb));
%                subject to
%                    x in X
%                    BarX.a' * x >= BarX.b;
%                    Bundle.const(1:Bundle.size) + Bundle.matrix(1:Bundle.size,:) * x <= LS * ones(Bundle.size, 1);
%            cvx_end
%            lb_t = cvx_optval;

        %specify the linear constrain matrix
        lb_prob.c = g;
        %specify the linear constrain matrix
        lb_prob.a = [-sparse(BarX.a'); sparse(Bundle.matrix(1:Bundle.size,:))];
        lb_prob.blc = [-inf; -inf*ones(Bundle.size,1)];
        lb_prob.buc = [-BarX.b; LS * ones(Bundle.size, 1) - Bundle.const(1:Bundle.size)];
        lb_prob.blx = zeros(domain.n,1); 
        %lb_prob.blx = -data.DomainType*ones(domain.n,1);
        lb_prob.bux = data.DomainType*ones(domain.n,1);
    else,

%            cvx_begin
%            cvx_quiet(true);
%                variables x(domain.n) t;
%                minimize t;
%                subject to
%                    x in X;
%                    BarX.a' * x >= BarX.b;
%                    Bundle.const(1:Bundle.size) + Bundle.matrix(1:Bundle.size,:) * x <= LS * ones(Bundle.size, 1);
%                    f + g' * (x - x_lb) < t;
%                    Bundle.const(1:Bundle.size) + Bundle.matrix(1:Bundle.size,:) * x < t;
%            cvx_end

        %specify the linear constrain matrix
        lb_prob.c = [zeros(1,domain.n),1]';

        %specify the linear constrain matrix
        lb_prob.a = [-sparse([BarX.a',0])]; 
        lb_prob.a = [lb_prob.a; sparse([Bundle.matrix(1:Bundle.size,:),zeros(Bundle.size,1)])];
        lb_prob.a = [lb_prob.a; [g',-1]];
        lb_prob.a = [lb_prob.a; sparse([Bundle.matrix(1:Bundle.size,:),-ones(Bundle.size,1)])];

        lb_prob.blc = [[]];
        lb_prob.buc = [-BarX.b; LS * ones(Bundle.size, 1) - Bundle.const(1:Bundle.size)];
        lb_prob.buc = [lb_prob.buc; -f + g' * x_lb; - Bundle.const(1:Bundle.size)];    
        
        lb_prob.blx = [zeros(domain.n,1); -inf];
        %lb_prob.blx = [-data.DomainType*ones(domain.n,1); -inf];
        lb_prob.bux = [data.DomainType*ones(domain.n,1); inf];
    end
end;

[r,lb_res] = mosekopt('minimize echo(0)', lb_prob );

%if strcmp(lb_res.sol.itr.solsta, 'DUAL_INFEASIBLE_CER'),
%    lb_t = -inf;
%else
%    lb_sol = lb_res.sol.itr.xx;
%    lb_mosek = f + g' * (lb_sol(1:domain.n) - x_lb);
%    lb_t = lb_mosek;
%end
    
lb_sol = lb_res.sol.itr.xx;
%lb_res.sol.itr.prosta
if r ~= 0 
%    if r == 4006 % MSK_RES_TRM_STALL
%       lb_sol = lb_res.sol.bas.xx;
%    else
%     disp(lb_res.rcodestr);
     lb_t = Inf; %%primal is infeasible
%    end 
else
    if strcmp(lb_res.sol.itr.prosta, 'PRIMAL_INFEASIBLE'),
        lb_t = Inf;
    else    
        lb_mosek = f + g' * (lb_sol(1:domain.n) - x_lb);
        lb_t = lb_mosek;
    end
end

