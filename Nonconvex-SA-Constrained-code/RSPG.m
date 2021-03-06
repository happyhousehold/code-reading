%% RSPG
m = data.m;
%N_iter = ceil(N_total/m);
N_total = data.m*N_iter;
fprintf('     \n');
fprintf('*** Running RSPG\n');
num_RSPG = num_RSPG+1;
if N_total <= 1000
    seed_val = data.seed;
elseif   N_total <= 5000
    seed_val = data.seed+20;
else
    seed_val = data.seed+40;
end
Run_number = 0;
j = 1;
temp = S;
S = 1;
Run_count = 0;
grad1_best  = zeros(1, Run_times);
grad2_best  = zeros(1, Run_times);
loss_best = zeros(1, Run_times);
obj_best = zeros(1, Run_times);
true_zero=zeros(1, Run_times);
false_zero=zeros(1, Run_times);
best_solution = zeros(Run_times, data.dim);
num_rnd = 1;
while j <= Run_times
    data.vali = ceil(N_iter/50);
    eval_xbar;
    j = j+1;
end
S = temp;
var_grad1 = var(grad1_best);
mean_grad1 = mean(grad1_best);
var_grad2 = var(grad2_best);
mean_grad2 = mean(grad2_best);
var_loss = var(loss_best);
mean_loss = mean(loss_best);
var_obj = var(obj_best);
mean_obj = mean(obj_best);
mean_true_zero = mean(true_zero);
mean_false_zero = mean(false_zero);
Results_showing;
best_solution_RSPG(count_RSPG+1:count_RSPG +Run_times,1) = N_iter;
best_solution_RSPG(count_RSPG+1:count_RSPG+Run_times,2:data.dim+1) = best_solution;
var_grad1_RSPG(num_RSPG,1) = N_iter;
var_grad1_RSPG(num_RSPG,2) = var_grad1;
var_grad2_RSPG(num_RSPG,1) = N_iter;
var_grad2_RSPG(num_RSPG,2) = var_grad2;
mean_grad1_RSPG(num_RSPG,1) = N_iter;
mean_grad1_RSPG(num_RSPG,2) = mean_grad1;
mean_grad2_RSPG(num_RSPG,1) = N_iter;
mean_grad2_RSPG(num_RSPG,2) = mean_grad2;
var_loss_RSPG(num_RSPG,1) = N_iter;
var_loss_RSPG(num_RSPG,2) = var_loss;
mean_loss_RSPG(num_RSPG,1) = N_iter;
mean_loss_RSPG(num_RSPG,2) = mean_loss;
var_obj_RSPG(num_RSPG,1) = N_iter;
var_obj_RSPG(num_RSPG,2) = var_obj;
mean_obj_RSPG(num_RSPG,1) = N_iter;
mean_obj_RSPG(num_RSPG,2) = mean_obj;
mean_true_zero_RSPG(num_RSPG,1) = N_iter;
mean_true_zero_RSPG(num_RSPG,2) = mean_true_zero;
mean_false_zero_RSPG(num_RSPG,1) = N_iter;
mean_false_zero_RSPG(num_RSPG,2) = mean_false_zero;
save(filename,'var_grad1_RSPG' ,'mean_grad1_RSPG','var_grad2_RSPG' ,'mean_grad2_RSPG','var_loss_RSPG','mean_loss_RSPG','var_obj_RSPG','mean_obj_RSPG','mean_true_zero_RSPG','mean_false_zero_RSPG','best_solution_RSPG', '-append')
frep=fopen(curr_path2,'a');
fprintf(frep,'     \n');
fprintf(frep,'*** Running RSPG\n');
Results_showing2;
fclose(frep);
count_RSPG = count_RSPG + Run_times;
