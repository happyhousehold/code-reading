%% Sample generation for SVM problem
data.dim = input('Enter the dimension of problem : ');
data.spr = input('Enter the sparsity of data set : ');
instance_matrix = ceil(sprand(N_initial,data.dim, data.spr));
z_sample = 2*rand(data.dim, 1)-ones(data.dim,1);
data.sample = z_sample;
label_vector = sign(instance_matrix*z_sample);
d = data.dim;
data.lambda = 0.01;
estimate_matrix = instance_matrix(1:N_initial,:);
estimate_lable = label_vector(1:N_initial);
z_ini = 2*rand(data.dim, 1)-ones(data.dim,1);
filename = ['SVM','-n',num2str(d),'.mat'];
data.seed = last_seed;
save (filename, 'data', 'z_ini' , 'z_sample', 'estimate_matrix', 'estimate_lable');
clear instance_matrix label_vector

