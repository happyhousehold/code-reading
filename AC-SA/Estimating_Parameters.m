load (filename);
N = N_initial;
z_sample = z_initial;
initial_solution = R_ini * z_ini;
A = (2/N)*x_initial'*x_initial + 2*lamda*eye(d);
% L is the lipschitz parameter
L = max(eig(A));
% mu is the strong convexity parameter
mu = min(eig(A));
