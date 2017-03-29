%% Sample generation for Trigonometric problem
data.dim = input('Enter the dimension of problem : ');
data.st = input('Enter the standard variation of noise : ');
z_ini = rand(data.dim,1);
n = data.dim;
filename = ['Trig','-n',num2str(n), '-st', num2str(st),'.mat'];
data.seed = last_seed;
save (filename, 'data', 'z_ini');