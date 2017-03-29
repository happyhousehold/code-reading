gamma = L_correction*min(1/(4*L*sqrt(data.dim+4)) , D_bar/sqrt(N_iter*sigma^2))/sqrt(data.dim+4);
t = 1;
all_iterates = zeros(r, data.dim);
while  t <= r
    zero_SGradient;
    z = z - gamma*GR;
    all_iterates(t,:) =z;
    t = t+1;
end