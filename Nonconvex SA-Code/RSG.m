%% RSG
fprintf('     \n');
fprintf('*** Running RSG\n');
num_RSG = num_RSG+1;
if N_iter <= 1000
    seed_val = data.seed;
elseif   N_iter <= 5000
    seed_val = data.seed+20;
else
    seed_val = data.seed+40;
end
Run_number = 0;
j = 1;
temp = S;
S = 1;
Run_count = 0;
grad_best  = zeros(1, Run_times);
loss_best = zeros(1, Run_times);
best_solution = zeros(Run_times, data.dim);
num_rnd = 1;
while j <= Run_times
    data.vali = ceil(N_iter/50);
    eval_xbar;
    j = j+1;
end
S = temp;
var_grad = var(grad_best);
mean_grad = mean(grad_best);
var_loss = var(loss_best);
mean_loss = mean(loss_best);
Results_showing;
best_solution_RSG(count_RSG+1:count_RSG +Run_times,1) = N_iter;
best_solution_RSG(count_RSG+1:count_RSG+Run_times,2:data.dim+1) = best_solution;
var_grad_RSG(num_RSG,1) = N_iter;
var_grad_RSG(num_RSG,2) = var_grad;
mean_grad_RSG(num_RSG,1) = N_iter;
mean_grad_RSG(num_RSG,2) = mean_grad;
var_loss_RSG(num_RSG,1) = N_iter;
var_loss_RSG(num_RSG,2) = var_loss;
mean_loss_RSG(num_RSG,1) = N_iter;
mean_loss_RSG(num_RSG,2) = mean_loss;
save(filename,'var_grad_RSG' ,'mean_grad_RSG','var_loss_RSG','mean_loss_RSG','best_solution_RSG', '-append')
frep=fopen(curr_path2,'a');
fprintf(frep,'     \n');
fprintf(frep,'*** Running RSG\n');
Results_showing2;
fclose(frep);
count_RSG = count_RSG + Run_times;
