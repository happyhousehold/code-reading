%%%% compute prox-mapping given by
%%%% s = \argmin_{x \in X} \langle p, x \rangle + \Omega(x)
%%%% where X is the standard spectahedron and 
%%%%%      \Omega(x) = Trace((x+\delta I) \ln (x+\deltaI)

function [s] = ProxmappingSimplex(data, domain, p)
n = domain.n;
R = domain.R;

mx = max(-p-1- (1e-10));
ds = exp(-p-1-(1e-10)-mx);
s = ds/sum(ds)*R;



