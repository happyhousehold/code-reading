%% Sample generation for Regression problem
fprintf('  \n');
d = input('Enter the dimension of problem : ');
data.dim = d;
data.st = input('Enter the standard variation of noise : ');
st = data.st;
data.spr = input('Enter the sparsity of data : ');
er_initial = normrnd(0,data.st,N_initial,1);
x_initial = sprand(N_initial, data.dim, data.spr);
z_sample = rand(data.dim, 1);
z_ini = rand(data.dim, 1);
y_initial = x_initial*z_sample + er_initial;
if strcmp(data.problem, 'Reg') == 1
    filename = ['Reg-n',num2str(d),'-st',num2str(st),'.mat'];
    data.lambda = 0.01;
else
    filename = ['Least_sqr-n',num2str(d),'-st',num2str(st),'.mat'];
    data.lambda = 0;
end
data.seed = last_seed;
save (filename, 'N_initial','data', 'er_initial', 'x_initial', 'z_sample', 'z_ini', 'y_initial');