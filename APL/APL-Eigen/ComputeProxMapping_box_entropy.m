function [sol] = ComputeProxMapping_box_entropy(domain, Bundle, data, BarX, f, g, x_lbt, c, LS)
%%%% use mosek conic optimier %%%%
clear map_prob;
clear map_res;
clear map_sol;
 
% the objective is given by
% \sum_i x_i / n log (x_i/n) - sum_i x_i/n(1+log c_i /n)
% set y = x / n, the objective becomes
% \sum_i y_i log(y_i) - sum_i y_i (1 + log c_i/n)

domega_c = 1 + log(c/domain.n);
map_prob.d = [sparse(domain.n,1);ones(domain.n,1)];
map_prob.c = [sparse(domain.n,1); -domega_c];
%map_prob.a = [ones(1,domain.n),sparse(1,domain.n)]; 
map_prob.a = [sparse(BarX.a'), sparse(1,domain.n)];
map_prob.a = [map_prob.a; g', sparse(1,domain.n)];
map_prob.a = [map_prob.a; sparse(Bundle.matrix(1:Bundle.size,:)), sparse(Bundle.size, domain.n)];
map_prob.a = [map_prob.a; sparse(eye(domain.n)), -domain.n*sparse(eye(domain.n))];
map_prob.blc = [BarX.b; -inf * ones(Bundle.size+1,1);sparse(domain.n,1)];
map_prob.buc = [inf; LS-f+g'*x_lbt; LS * ones(Bundle.size, 1)-Bundle.const(1:Bundle.size);sparse(domain.n,1)];

map_prob.blx = [zeros(domain.n,1);-inf * ones(domain.n+1,1)]; 
map_prob.bux = [ones(domain.n,1);inf * ones(domain.n+1,1)];

map_res = mskenopt(map_prob.d, map_prob.c, map_prob.a,map_prob.blc,map_prob.buc,[],'minimize echo(0)');
%map_res = mskenopt(map_prob.d, map_prob.c, map_prob.a,map_prob.blc,map_prob.buc);

map_sol = map_res.sol.itr.xx;
    
sol = map_sol(1:domain.n);