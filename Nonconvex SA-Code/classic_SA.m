t = 1;
D_bar = sqrt(2*e_0/L);
gamma = L_correction* min(1/L ,D_bar/sqrt(N_iter*sigma^2));
all_iterates = zeros(r, data.dim);
while  t <= r
    SGradient;
    z = z - gamma*GR;
    all_iterates(t,:) = z;
    t = t+1;
end