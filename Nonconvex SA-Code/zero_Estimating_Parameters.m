e_0 = 130;
data.mu = 0.0025;
L1 = zeros(20,1);
for i = 1 : 20
    x(1)= randi(100);
    x(2) = x(1) + randi(100);
    y(1)= randi(100);
    y(2) = y(1) + randi(100);
    stch_grad1 = zeros(20,2);
    stch_grad2 = zeros(20,2);
    L11 = zeros(20,2);
    for j = 1 : 20
        z = x';
        zero_SGradient;
        G1 = GR;
        z = y';
        zero_SGradient;
        G2 = GR;
        stch_grad1(2*j-1,:) = G1;
        stch_grad2(2*j,:) = G2;
        L11(j,:) = G1 - G2;
    end
    sig(2*i-1) = sqrt(sum(var(stch_grad1)));
    sig(2*i) = sqrt(sum(var(stch_grad2)));
    L1(i) = norm(mean(L11))/norm(x-y);
end
L = max(L1);
sigma = max(sig);
D_bar = 2*e_0/L;