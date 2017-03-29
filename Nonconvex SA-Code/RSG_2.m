%% 2-RSG
fprintf('     \n');
fprintf('*** Running 2-RSG\n');
num_2RSG = num_2RSG+1;
if N_iter <= 200
    seed_val = data.seed+60;
elseif   N_iter <= 1000
    seed_val = data.seed+160;
else
    seed_val = data.seed+260;
end
Run_number = 0;
j = 1;
grad_best  = zeros(1, Run_times);
loss_best  = zeros(1, Run_times);
best_solution = zeros(Run_times, data.dim);
num_rnd = 1;
while j <= Run_times
    data.vali = N_iter/2;
    eval_xbar;
    j = j+1;
end
var_grad = var(grad_best);
mean_grad = mean(grad_best);
var_loss  = var(loss_best);
mean_loss = mean(loss_best);
Results_showing;
best_solution_2RSG(count_2RSG+1:count_2RSG +Run_times,1) = N_iter;
best_solution_2RSG(count_2RSG+1:count_2RSG+Run_times,2:data.dim+1) = best_solution;
var_grad_2RSG(num_2RSG,1) = N_iter;
var_grad_2RSG(num_2RSG,2) = var_grad;
mean_grad_2RSG(num_2RSG,1) = N_iter;
mean_grad_2RSG(num_2RSG,2) = mean_grad;
var_loss_2RSG(num_2RSG,1) = N_iter;
var_loss_2RSG(num_2RSG,2) = var_loss;
mean_loss_2RSG(num_2RSG,1) = N_iter;
mean_loss_2RSG(num_2RSG,2) = mean_loss;
save(filename,'var_grad_2RSG' ,'mean_grad_2RSG','var_loss_2RSG','mean_loss_2RSG','best_solution_2RSG', '-append')
frep=fopen(curr_path2,'a');
fprintf(frep,'     \n');
fprintf(frep,'*** Running 2-RSG\n');
Results_showing2;
fclose(frep);
count_2RSG = count_2RSG + Run_times;
