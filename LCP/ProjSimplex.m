%%%% compute prox-mapping given by
%%%% s = \argmin_{x \in X} \langle p, x \rangle + \|x\|^2/2
%%%% where X is the standard spectahedron 

function [s] = ProjSimplex(data, domain, p)
n = domain.n;
R = domain.R;
%%% do eigenvalue decomposition for p.
%[V, xi] = eig(full(reshape(p, n,n)));
%q = diag(xi);

% now we need to solve the problem of
% argmin_{x in Simplex} \langle q, x \rangle + \|x\|^2/2

%nq = sort(-p);

%for i = 1: n-1
%    lambda = (sum(nq(i:n))-R)/(n-i+1);
%    if i==1 
%        if lambda <= nq(i)
%            break;
%        end
%    else
%        if lambda > nq(i-1) & lambda <= nq(i)
%            break
%        end
%    end
    
   
%end

%s = max(0, nq - lambda);
x = zeros(n,1);
sumx = 0;
for i = 1:n 
   sumx = sumx + p(i);
end
%        // find a lambda 
%        // sort those xi's
sorted = sort(p);

for i = 1:n
   lambda =  -(R + sumx) / (n - i +1);
   if lambda < - sorted(n - i +1)
        break;
   end
   sumx = sumx - sorted(n - i +1);
end
        
x= R* max(-p-lambda, 0);

      
s=x;
