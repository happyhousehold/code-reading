function [f,g,f0] = Oracle(data, x)
%%%%%%%%
% Let p = [zeros(data.n2-data.dim_xi,1);ones(dim_xi,1)]
% the primal form of the second stage problem is given by:
% c'x + min p'f
% s.t  
%      A_1 f <= -A_2 x + b
%          f >= 0
% Hence the dual is given by
% v(x,\xi):= max (-A_2 x + b)' y
% s.t. 
%      A_1' y ?? p
%           y1 ?? 0
% The sample average of the second stage problem is given by
%
% F(x) = c'x + 1/n \sum_{i=1}^ns V(x,xi_i),
%
% and the gradient is given by
% c - 1/n A_2' * \sum_{i=1}^ns y2(\xi_i)
%
% One can also smooth V(x) by defining
%F_\mu(x) = c'x + 1/n \sum_i  V_\mu(x,\xi_i)
%where V_\mu(x,xi) = max (-A_2 x + b)'y - \mu \|y\|^2 /2
%                    s.t.   A_1' y ?? p, y1 free, y2 <= 0



f = -inf;
sumf = 0;

g = zeros(data.n,1);
sumg = g;

f0 = -inf;
sumf0 = 0;

for i = 1:data.ns, %% totally ns senarios
    %% compute f0
    data.org_dual.c = data.rand_second_rhs(i,:)' - data.A2(:, 1:data.n1) * x;
    % solve the orginal dual
    [r,od_res] = mosekopt('maximize echo(0)', data.org_dual);
    if strcmp(od_res.sol.itr.prosta, 'DUAL_INFEASIBLE'),
        sumf0 = inf;
        break;
    else        
        od_sol = od_res.sol.itr.xx;
        sumf0 = sumf0 + data.org_dual.c'*od_sol;
    end
    
    if data.mu == 0 % no smooth
        sumg = sumg - data.A2(:, 1:data.n1)' * od_sol;
    else % solve the smoothed dual
        data.sm_dual.q = data.mu * eye(data.m2) /2;
        data.sm_dual.c = -data.org_dual.c;
        [r,sd_res] = mosekopt('minimize echo(0)', data.sm_dual);
        sd_sol = sd_res.sol.itr.xx;
        sumf = sumf + data.org_dual.c'*sd_sol - data.mu / 2 * sd_sol' * sd_sol;
        
        sumg = sumg  - data.A2(:, 1:data.n1)' * sd_sol;
    end
end

f0 = data.master.c'*x + sumf0 / data.ns;
g = data.master.c + sumg / data.ns;
if data.mu == 0,
    f = f0;
else 
    f = data.master.c'*x + sumf / data.ns;
end

