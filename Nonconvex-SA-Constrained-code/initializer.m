curr_path = pwd;
curr_path = [curr_path,'\Result'];
curr_path2 = [curr_path,'\Results.txt'];
frep = fopen(curr_path2,'a');
fprintf(frep,'     \n');
filename = ' ';
running = 0;
S = 5;
N_iter = 1000;
N_initial = 200;
N_training = N_initial;
%lambda = sqrt(2*log2(data.dim));
lambda=0.01;
N_vali = 50000;
Run_times = 20;
R_ini = 1;
L_correction = 1;
data.vali = N_vali;
%data.lambda = lambda;
data.lambda1 = 1;
data.lambda2 = 0.5;
data.lambda3 = 0.5;
parallel_num = 1;
ep = 0.1;
p = 0.5;
data.a=3.7;
N_total=25000;
data.threshold = 0.02;
alpha=1;
Run_number = 0;
count_RSPG = 0;
num_RSPG = 0;
count_2RSPG = 0;
num_2RSPG = 0;
count_2RSPG_V = 0;
num_2RSPG_V = 0;
count_RSG = 0;
num_RSG = 0;
count_2RSG = 0;
num_2RSG = 0;
count_2RSG_V = 0;
num_2RSG_V = 0;
count_MD_SA = 0;
num_MD_SA = 0;
seed_val = 0;
num_rnd =1;
generating_rnd_num;
