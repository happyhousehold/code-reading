%% Estimating problem parameters
initial_solution = R_ini * z_ini;
% data.lambda  = sqrt(2*log2(data.dim));
data.lambda = 0.01;
%% Non-convex S3VM problem
sig = zeros(1,5);
A = zeros(data.dim);
for i =1 : N_initial
    A = A + 2*data.lambda2*estimate_lable1(i)*repmat(estimate_matrix1(i,:),data.dim,1).*repmat(estimate_matrix1(i,:)',1,data.dim)+20*data.lambda3*repmat(estimate_matrix2(i,:),data.dim,1).*repmat(estimate_matrix2(i,:)',1,data.dim);
end
A = A/N_initial+data.lambda1*eye(data.dim);
L = max(real(eig(A)));
sig = zeros(1, 10);
for i = 1 : 5
    initial = 2*rand(data.dim, 1) - rand(data.dim, 1);
    stch_grad = zeros(20,data.dim);
    for j = 1 : 20
        GR_b = 1-estimate_lable1((i-1)*20+j)*(estimate_matrix1((i-1)*20+j,:)*initial);
        GR_b = 2*GR_b*(GR_b>0);
        GR2 = [-estimate_lable1((i-1)*20+j)*GR_b ; (estimate_matrix1((i-1)*20+j,2:end)*(-estimate_lable1((i-1)*20+j)*GR_b))'];
        GR_b = -10*(estimate_matrix2((i-1)*20+j,:)*initial)*exp(-5*(estimate_matrix2((i-1)*20+j,:)*initial).^2);
        GR3 = [GR_b ; (estimate_matrix2((i-1)*20+j,2:end)*GR_b)'];
        stch_grad(j,:)= (data.lambda2*GR2+ data.lambda3*GR3)';
    end
    sig(i) = sqrt(sum(var(stch_grad)));
end
sigma = max(sig);
data.vali = 5000;
data.gamma  =  alpha/(2*L);
[loss_0, e_0 , g_0] = evaluation_oracle(initial_solution,data);
data.m = ceil(max(1, sqrt(2.5*N_total)*sigma/(2*sqrt(2*L*e_0))));
data.ratio  = sum(estimate_lable1>0)/N_initial;
data.toler = 0.1;

data.vali=75000;
