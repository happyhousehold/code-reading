%%%% compute prox-mapping given by
%%%% s = \argmin_{x \in X} \langle p, x \rangle + \|x\|^2/2
%%%% where X is the standard spectahedron 

function s = ProjSpectahedron(data, domain, p)
n = domain.n;
R = domain.R;
%%% do eigenvalue decomposition for p.
sp = reshape(p, n,n);
sp = (sp + sp')/2;
[V, xi] = eig(full(sp));
q = diag(xi);

% now we need to solve the problem of
% argmin_{x in Simplex} \langle q, x \rangle + \|x\|^2/2

%nq = sort(-q);

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

%dx = max(0, nq - lambda);

x = zeros(n,1);
sumx = 0;
for i = 1:n 
   sumx = sumx + q(i);
end
%        // find a lambda 
%        // sort those xi's
sorted = sort(q);

for i = 1:n
   lambda = -(R + sumx) / (n - i +1);
   if lambda < - sorted(n - i +1)
        break;
   end
   sumx = sumx - sorted(n - i+1);
end
        
%x= R* max(-q-lambda, 0);
x = -q-lambda;
for i = 1:n
    if x(i) < 0
        x(i) = 0;
    end
end

%if min(eig(V*diag(x)*V')) < -1e-8
%    disp(sprintf('error with projection1'));
    
%    if min(x) < -1e-8
%        disp(sprintf('error with projection0'));
%    end
%end
s  = reshape(V * diag(x) * V', n^2, 1);

%if min(eig(reshape(s,n,n))) < -1e-8
%    disp(sprintf('error with projection2'));
%end
%s  = reshape(V' * diag(x) * V, n^2, 1);


