function [sol] = ComputeProxMapping_box(domain, Bundle, data, BarX, f, g, x_lbt, c, LS)
%%%% use mosek conic optimier %%%%
clear map_prob;
clear map_res;
clear map_sol;
 
%        cvx_begin
%        cvx_quiet(true);
%            variables x(domain.n), y(domain.n),t1;
%            minimize t2;
%            subject to
%                x - y  = c;
%                BarX.a' * x >= BarX.b;
%                Bundle.const(1:Bundle.size) + Bundle.matrix(1:Bundle.size, :) * x <= LS * ones(Bundle.size, 1);
%                norm(y) <= t1;
%                x in box;
%        cvx_end

    map_prob.c = [zeros(1, 2*domain.n),1]'; % first n entries are x, second n entries are y, last two entries: t1,t2
    map_prob.a = [sparse(eye(domain.n)),-sparse(eye(domain.n)),sparse(domain.n,1)];
    map_prob.a = [map_prob.a; sparse(BarX.a'), sparse(1,domain.n+1)];
    map_prob.a = [map_prob.a; sparse(Bundle.matrix(1:Bundle.size, :)), sparse(Bundle.size,domain.n+1)];
    map_prob.blc = [c;BarX.b;-inf*ones(Bundle.size,1)];
    map_prob.buc = [c;inf; LS * ones(Bundle.size, 1)-Bundle.const(1:Bundle.size)];
    map_prob.blx = [zeros(domain.n,1);-inf * ones(domain.n+1,1)]; 
  %  map_prob.blx = [-data.DomainType*ones(domain.n,1);-inf * ones(domain.n+1,1)];
    map_prob.bux = [data.DomainType*ones(domain.n,1);inf * ones(domain.n+1,1)];

    map_prob.cones = cell(1,1);
      
    map_prob.cones{1}.type = 'MSK_CT_QUAD';
    map_prob.cones{1}.sub  = [2*domain.n+1, [domain.n+1:1:2*domain.n]];

[r,map_res] = mosekopt('minimize echo(0)', map_prob );

map_sol = map_res.sol.itr.xx;
    
sol = map_sol(1:domain.n);