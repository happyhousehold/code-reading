t = 1;
D_bar = sqrt(2*e_0/L);
if m == 1
    gamma = min(1/L ,D_bar/sqrt(N_iter*sigma^2));
else
    gamma = alpha/(2*L);
end
all_iterates = zeros(r, data.dim);
all_iterates(1,:) = z;
while  t < r
    SGradient;
    z = z - gamma*GR;
    %S3VM
    if z(1) > 2*data.ratio-1+data.toler
        z(1) = 2*data.ratio-1+data.toler;
    elseif z(1) < 2*data.ratio-1-data.toler
        z(1) = 2*data.ratio-1-data.toler;
    end    
        
    %Lasso
    %z  =  sign(z - gamma*GR).*max(abs(z - gamma*GR)-gamma*data.lambda,0);
    all_iterates(t+1,:) = z;
    t = t+1;
end