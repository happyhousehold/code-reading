%% RSGF
fprintf('     \n');
fprintf('*** Running RSGF\n');
num_RSGF = num_RSGF+1;
if N_iter <= 1000
    seed_val = data.seed;
else
    seed_val = data.seed+10;
end
Run_number = 0;
total_run = 0;
j = 1;
temp = S;
S = 1;
Run_count = 0;
grad_best  = zeros(1, Run_times);
loss_best = zeros(1, Run_times);
best_solution = zeros(Run_times, data.dim);
num_rnd = 1;
while j <= Run_times
    data.vali = ceil(N_iter/100);
    zero_eval_xbar;
    j = j+1;
end
S = temp;
var_grad = var(grad_best);
mean_grad = mean(grad_best);
var_loss = var(loss_best);
mean_loss = mean(loss_best);
zero_Results_showing;
best_solution_RSGF(count_RSGF+1:count_RSGF +Run_times,1) = N_iter;
best_solution_RSGF(count_RSGF+1:count_RSGF+Run_times,2:data.dim+1) = best_solution;
var_grad_RSGF(num_RSGF,1) = N_iter;
var_grad_RSGF(num_RSGF,2) = var_grad;
mean_grad_RSGF(num_RSGF,1) = N_iter;
mean_grad_RSGF(num_RSGF,2) = mean_grad;
var_loss_RSGF(num_RSGF,1) = N_iter;
var_loss_RSGF(num_RSGF,2) = var_loss;
mean_loss_RSGF(num_RSGF,1) = N_iter;
mean_loss_RSGF(num_RSGF,2) = mean_loss;
save(filename, 'initial_solution','var_grad_RSGF' ,'mean_grad_RSGF','var_loss_RSGF','mean_loss_RSGF','best_solution_RSGF', '-append')
frep=fopen(curr_path2,'a');
fprintf(frep,'     \n');
fprintf(frep,'*** Running RSGF\n');
zero_Results_showing2;
fclose(frep);
count_RSGF = count_RSGF + Run_times;
