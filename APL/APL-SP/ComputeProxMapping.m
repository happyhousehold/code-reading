function [sol] = ComputeProxMapping(domain, Bundle, data, BarX, f, g, x_lbt, c, LS)
%%%% use mosek conic optimier %%%%
clear map_prob;
clear map_res;
clear map_sol;

%%% numerical statbility problems may exist if we change the toteranc
%%% in the experiment, we use the defalut settting for the instance
%%% 20-strom  (100) but use the following setting for 20-term (50)

% Primal feasibility tolerance for the primal solution
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1.0e-20;
% Dual feasibility tolerance for the dual solution
param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = 1.0e-20;
% Relative primal - dual gap tolerance .
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0e-20;
%        cvx_begin
%        cvx_quiet(true);
%            variables x(domain.n), y(domain.n),t1;
%            minimize t2;
%            subject to
%                x \in X
%                x - y  = c;
%                BarX.a' * x >= BarX.b;
%                Bundle.const(1:Bundle.size) + Bundle.matrix(1:Bundle.size, :) * x <= LS * ones(Bundle.size, 1);
%                norm(y) <= t1;
%        cvx_end

    map_prob.c = [zeros(1, 2*domain.n),1]'; % first n entries are x, second n entries are y, last two entries: t1,t2
    map_prob.a = [data.master.a, sparse(data.m1, domain.n+1)];
    map_prob.a = [map_prob.a; sparse(eye(domain.n)),-sparse(eye(domain.n)),sparse(domain.n,1)];
    map_prob.a = [map_prob.a; sparse(BarX.a'), sparse(1,domain.n+1)];
    map_prob.a = [map_prob.a; sparse(Bundle.matrix(1:Bundle.size, :)), sparse(Bundle.size,domain.n+1)];
    map_prob.blc = [data.master.blc; c; BarX.b;-inf*ones(Bundle.size,1)];
    map_prob.buc = [data.master.buc; c; inf; LS * ones(Bundle.size, 1)-Bundle.const(1:Bundle.size)];
    map_prob.blx = [data.master.blx;-inf * ones(domain.n+1,1)]; 
  %  map_prob.blx = [-data.DomainType*ones(domain.n,1);-inf * ones(domain.n+1,1)];
    map_prob.bux = [data.master.bux;inf * ones(domain.n+1,1)];

    map_prob.cones = cell(1,1);
      
    map_prob.cones{1}.type = 'MSK_CT_QUAD';
    map_prob.cones{1}.sub  = [2*domain.n+1, [domain.n+1:1:2*domain.n]];

%[r,map_res] = mosekopt('minimize echo(0)', map_prob, param );
[r,map_res] = mosekopt('minimize echo(0)', map_prob );

map_sol = map_res.sol.itr.xx;

sol = map_sol(1:domain.n);