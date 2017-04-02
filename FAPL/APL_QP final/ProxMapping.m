function [sol] = ProxMapping(data, Bundle, control, BarX, c, LS)
%%%% use mosek conic optimier %%%%
clear map_prob;
clear map_res;
clear map_sol;
         
% How to change the parameters that controls
% the accuracy of a solution computed by the conic
% optimizer .
param = [];
% Primal feasibility tolerance for the primal solution
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1.0e-20;
% Dual feasibility tolerance for the dual solution
param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = 1.0e-20;
% Relative primal - dual gap tolerance .
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0e-20;
domain.n= data.n;
domain.R= data.R;
data.mu=0;
if data.mu < 1e-2 
    
%        cvx_begin
%        cvx_quiet(true);
%            variables x(domain.n), y(domain.n),t1,t2;
%            minimize t2;
%            subject to
%                t1 = domain.R;
%                x - y  = c;
%                BarX.a' * x >= BarX.b;
%                Bundle.const(1:Bundle.size) + Bundle.matrix(1:Bundle.size, :) * x <= LS * ones(Bundle.size, 1);
%                norm(x) <= t1;
%                norm(y) <= t2;
%        cvx_end

    map_prob.c = [zeros(1, 2*domain.n),0,1]'; % first n entries are x, second n entries are y, last two entries: t1,t2
    map_prob.a = [sparse(1,2*domain.n),sparse([1,0])];
    map_prob.a = [map_prob.a; sparse(eye(domain.n)),-sparse(eye(domain.n)),sparse(domain.n,2)];
    map_prob.a = [map_prob.a; sparse(BarX.a'), sparse(1,domain.n+2)];
    map_prob.a = [map_prob.a; sparse(Bundle.matrix(1:Bundle.size, :)), sparse(Bundle.size,domain.n+2)];
    map_prob.blc = [domain.R;c;BarX.b;-inf*ones(Bundle.size,1)];
    map_prob.buc = [domain.R;c;inf; LS * ones(Bundle.size, 1)-Bundle.const(1:Bundle.size)];
    
    map_prob.cones = cell(2,1);
    map_prob.cones{1}.type = 'MSK_CT_QUAD';
    map_prob.cones{1}.sub  = [2*domain.n+1, [1:1:domain.n]];
    
    map_prob.cones{2}.type = 'MSK_CT_QUAD';
    map_prob.cones{2}.sub  = [2*domain.n+2, [domain.n+1:1:2*domain.n]];
else %% strongly convex case  
%        cvx_begin
%        cvx_quiet(true);
%            variables x(domain.n), y(domain.n), z(domain.n), t1,t2,t3;
%            minimize t2;
%            subject to
%                t1 = domain.R;
%                t3 = BarX.newR;
%                x - y  = c;
%                x - z = BarX.center
%                BarX.a' * x >= BarX.b;
%                Bundle.const(1:Bundle.size) + Bundle.matrix(1:Bundle.size, :) * x <= LS * ones(Bundle.size, 1);
%                norm(x) <= t1;
%                norm(y) <= t2;
%                norm(z) <= t3;
%        cvx_end

    map_prob.c = [zeros(1, 3*domain.n),0,1,0]'; 
    map_prob.a = [sparse(1,3*domain.n),sparse([1,0,0])];
    map_prob.a = [map_prob.a; sparse([0,0,1])];
    map_prob.a = [map_prob.a; sparse(eye(domain.n)),-sparse(eye(domain.n)),sparse(domain.n,domain.n+3)];
    map_prob.a = [map_prob.a; sparse(eye(domain.n)),sparse(domain.n,domain.n), -sparse(eye(domain.n)),sparse(domain.n,3)];
    
    map_prob.a = [map_prob.a; sparse(BarX.a'), sparse(1,2*domain.n+3)];
    map_prob.a = [map_prob.a; sparse(Bundle.matrix(1:Bundle.size, :)), sparse(Bundle.size,2*domain.n+3)];
    map_prob.blc = [domain.R;BarX.newR;c;BarX.center;BarX.b;-inf*ones(Bundle.size,1)];
    map_prob.buc = [domain.R;BarX.newR;c;BarX.center;inf; LS * ones(Bundle.size, 1)-Bundle.const(1:Bundle.size)];
    
    map_prob.cones = cell(3,1);
    map_prob.cones{1}.type = 'MSK_CT_QUAD';
    map_prob.cones{1}.sub  = [3*domain.n+1, [1:1:domain.n]];
    
    map_prob.cones{2}.type = 'MSK_CT_QUAD';
    map_prob.cones{2}.sub  = [3*domain.n+2, [domain.n+1:1:2*domain.n]];
    
    map_prob.cones{3}.type = 'MSK_CT_QUAD';
    map_prob.cones{3}.sub  = [3*domain.n+3, [2*domain.n+1:1:3*domain.n]];
end;

%[r,map_res] = mosekopt('minimize echo(0)', map_prob, param );
[r,map_res] = mosekopt('minimize echo(0)', map_prob );
map_sol = map_res.sol.itr.xx;
    
%if r ~= 0 
%    if r == 4006 % MSK_RES_TRM_STALL
%        map_sol = map_res.sol.bas.xx;
%    else
%        disp(map_res.rcodestr);
%    end
%end
sol = map_sol(1:domain.n);