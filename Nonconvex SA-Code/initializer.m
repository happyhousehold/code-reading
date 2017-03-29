filename = ' ';
running = 0;
S = 5;
N_iter = 200;
N_initial = 200;
N_training = N_initial;
if (strcmp(data.problem, 'Least_sqr') == 1) || (strcmp(data.problem, 'SVM') == 1)
    lambda = 0.01;
else
    lambda = 0;
end
N_vali = 75000;
Run_times = 20;
if (strcmp(data.problem, 'Least_sqr') == 1) || (strcmp(data.problem, 'Reg') == 1)
    R_ini = 0;
elseif (strcmp(data.problem, 'SVM') == 1)
    R_ini = 5;
elseif (strcmp(data.problem, 'Trig') == 1)
    R_ini = 1;
end
L_correction = 1;
data.vali = N_vali;
data.lambda = lambda;
parallel_num = 1;
ep = 0.1;
p = 0.5;
Run_number = 0;
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