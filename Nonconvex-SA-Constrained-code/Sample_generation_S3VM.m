%% Sample generation for SVM problem
data.dim = input('Enter the dimension of problem : ');
data.spr = input('Enter the sparsity of data set : ');
%instance_matrix = ceil(sprand(N_initial,data.dim, data.spr));
instance_matrix  = [ones(2*N_initial,1),sprandn(2*N_initial,data.dim-1,data.spr)];
%z_sample = sprand(data.dim, 1, data.spr)-sprand(data.dim, 1, data.spr);
z_sample = 2*rand(data.dim, 1)-ones(data.dim,1);
data.sample = z_sample;
label_vector = sign(instance_matrix(1:N_initial,:)*z_sample);
d = data.dim;
%sz = data.size;
estimate_matrix1 = instance_matrix(1:N_initial,:);
estimate_lable1 = label_vector;
estimate_matrix2 = instance_matrix(N_initial+1:end,:);
z_ini = 2*rand(data.dim, 1)-ones(data.dim,1);
data.test_matrix = [ones(2*N_vali,1),sprandn(2*N_vali,data.dim-1,data.spr)];
data.test_lable = sign(data.test_matrix(1:N_vali,:)*data.sample);
filename = ['S3VM','-n',num2str(d),'.mat'];
data.seed = last_seed;
save (filename, 'data', 'z_ini' , 'z_sample', 'estimate_matrix1','estimate_matrix2', 'estimate_lable1');
clear instance_matrix label_vector

