function [xt, lb] = LMO(gyt, domain)

%%%% determine the number of nonzeros in the solution
%%%% dependent on the domain type, generate an optimal solution
type = domain.type;
R = domain.R;
n = domain.n;

if type == 0 %% standard box intersected with simplex
    xt = zeros(n,1);
    [S_gyt, I_gyt] = sort(gyt);
    sum = 0;
    for i=1:n,
        if S_gyt(i) >= 0, %% the remaining entries must be >=0
            break;
        end
        
        if sum + 1 <= R
            xt(I_gyt(i)) = 1;
            sum = sum+1;
        else
            xt(I_gyt(i)) = R - sum;
            sum = R;
        end
        
        if sum == R,
            break;
        end
    end
    lb = gyt' * xt;
    
    
    %%% just for the purpose of testing
 %  xt1 = zeros(domain.n,1);
 %  xt1 = (gyt < 0); 
 %  lb1 = gyt'* xt1;
   
 %  if lb ~= lb1
 %      disp(sprintf('error,%.6e', norm(xt-xt1)));
 %  end
   
end

if type == 1 %% box
   xt = zeros(domain.n,1);
   xt = domain.R *(gyt < 0); 
   lb = gyt'* xt;
%   xt = domain.R *(gyt < 0)- domain.R*(gyt >= 0); 
%   lb = gyt'* xt;
end

if type == 2 %% simplex
   xt = zeros(domain.n,1);
   [v,in] = min(gyt); 
   xt(in) = R;
   lb = v*R;
end

if type == 3 %% ball
    xt = zeros(domain.n,1);
    v = norm(gyt);
    lb = -R* v;
    xt = -R * gyt / v
end
    
if type == 4 %% spectahedron
    %% find the minimum eigenvector/value for gyt
    opts.issym = 1;
    opts.disp = 0;
    [V,D] = eigs(sparse(reshape(gyt,n,n)), 1, 'sa', opts);
 %  [V,D] = myEig(reshape(-gyt,n,n), n, 1e-6);
 %  V = - V;
 %  D = - D;
    xt = reshape(R*V*V', n^2, 1);
    lb = D*R;
end

