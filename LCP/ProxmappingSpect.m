%%%% compute prox-mapping given by
%%%% s = \argmin_{x \in X} \langle p, x \rangle + \Omega(x)
%%%% where X is the standard spectahedron and 
%%%%%      \Omega(x) = Trace((x+\delta I) \ln (x+\deltaI)

function [s] = ProxmappingSpect(data, domain, p)
n = domain.n;
R = domain.R;
%%% do eigenvalue decomposition for p.
[V, xi] = eig(full(reshape(p, n,n)));

mx = max(diag(-xi-1- (1e-10)));
ds = exp(-xi-1-(1e-10)-mx);
ds = ds/sum(ds)*R;

s  = reshape(V * diag(ds) * V', n^2, 1);


