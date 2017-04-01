function data = NewLinearSystem(n, density, cn, R)

data.A = sprandsym(n,density,1.0/cn,1);
%data.A = rand(n,n);
tx = rand(n,1);
data.sol = R * tx / norm(tx);
data.b = data.A * data.sol;
v = eig(data.A' * data.A);
data.L = max(v);
data.mu = min(v);
