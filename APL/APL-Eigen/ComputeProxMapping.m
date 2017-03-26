function [sol] = ComputeProxMapping(domain, Bundle, data, BarX, f, g, x_lbt, c, LS)
%%%% use mosek conic optimier %%%%
clear map_prob;
clear map_res;
clear map_sol;

%Bnd = domain.n;
Bnd = 1;

% cvx_begin
%     cvx_quiet(true);
%         variable x(domain.n);
%         minimize -sum(entr(x+delta) -domega_c.*x);
%         subject to
%             x == simplex(domain.n);
%             BarX.a' * x >= BarX.b;
% %              f + g' * (x - x_lbt) <= LS;
%             Bundle.const(1:Bundle.size) + Bundle.matrix(1:Bundle.size, :) * x <= LS * ones(Bundle.size, 1);
%     cvx_end
x = zeros(domain.n,1);
%delta = 1.0e-16 / domain.n;
delta = 0;
domega_c = 1 + log(c + delta);
map_prob.d = ones(domain.n,1);
map_prob.c = -domega_c;
map_prob.a = [ones(1,domain.n); sparse(BarX.a'); g'; sparse(Bundle.matrix(1:Bundle.size,:))];
map_prob.blc = [Bnd; BarX.b; -inf * ones(Bundle.size+1,1)];
map_prob.buc = [Bnd; inf; LS-f+g'*x_lbt; LS * ones(Bundle.size, 1)-Bundle.const(1:Bundle.size)];

map_res = mskenopt(map_prob.d, map_prob.c, map_prob.a,map_prob.blc,map_prob.buc,[],'minimize echo(0)');
%map_res = mskenopt(map_prob.d, map_prob.c, map_prob.a,map_prob.blc,map_prob.buc);

map_sol = map_res.sol.itr.xx;
    
sol = map_sol(1:domain.n);