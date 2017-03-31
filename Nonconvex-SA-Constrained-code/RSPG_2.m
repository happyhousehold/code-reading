%% 2-RSPG
m = data.m;
%N_iter = ceil(N_total/m);
fprintf('     \n');
fprintf('*** Running 2-RSPG\n');
num_2RSPG = num_2RSPG+1;
if N_total <= 200
    seed_val = data.seed+60;
elseif   N_total <= 1000
    seed_val = data.seed+160;
else
    seed_val = data.seed+260;
end
Run_number = 0;
j = 1;
grad1_best  = zeros(1, Run_times);
grad2_best  = zeros(1, Run_times);
loss_best  = zeros(1, Run_times);
obj_best  = zeros(1, Run_times);
true_zero=zeros(1, Run_times);
false_zero=zeros(1, Run_times);
best_solution = zeros(Run_times, data.dim);
num_rnd = 1;
while j <= Run_times
    data.vali = ceil(N_iter/2);
    eval_xbar;
    j = j+1;
end
var_grad1 = var(grad1_best);
mean_grad1 = mean(grad1_best);
var_grad2 = var(grad2_best);
mean_grad2 = mean(grad2_best);
var_loss  = var(loss_best);
mean_loss = mean(loss_best);
var_obj  = var(obj_best);
mean_obj = mean(obj_best);
mean_true_zero = mean(true_zero);
mean_false_zero = mean(false_zero);
Results_showing;
best_solution_2RSPG(count_2RSPG+1:count_2RSPG +Run_times,1) = N_iter;
best_solution_2RSPG(count_2RSPG+1:count_2RSPG+Run_times,2:data.dim+1) = best_solution;
var_grad1_2RSPG(num_2RSPG,1) = N_iter;
var_grad1_2RSPG(num_2RSPG,2) = var_grad1;
var_grad2_2RSPG(num_2RSPG,1) = N_iter;
var_grad2_2RSPG(num_2RSPG,2) = var_grad2;
mean_grad1_2RSPG(num_2RSPG,1) = N_iter;
mean_grad1_2RSPG(num_2RSPG,2) = mean_grad1;
mean_grad2_2RSPG(num_2RSPG,1) = N_iter;
mean_grad2_2RSPG(num_2RSPG,2) = mean_grad2;
var_loss_2RSPG(num_2RSPG,1) = N_iter;
var_loss_2RSPG(num_2RSPG,2) = var_loss;
mean_loss_2RSPG(num_2RSPG,1) = N_iter;
mean_loss_2RSPG(num_2RSPG,2) = mean_loss;
var_obj_2RSPG(num_2RSPG,1) = N_iter;
var_obj_2RSPG(num_2RSPG,2) = var_obj;
mean_obj_2RSPG(num_2RSPG,1) = N_iter;
mean_obj_2RSPG(num_2RSPG,2) = mean_obj;
mean_true_zero_2RSPG(num_2RSPG,1) = N_iter;
mean_true_zero_2RSPG(num_2RSPG,2) = mean_true_zero;
mean_false_zero_2RSPG(num_2RSPG,1) = N_iter;
mean_false_zero_2RSPG(num_2RSPG,2) = mean_false_zero;
save(filename,'var_grad1_2RSPG' ,'mean_grad1_2RSPG','var_grad2_2RSPG' ,'mean_grad2_2RSPG','var_loss_2RSPG','mean_loss_2RSPG','var_obj_2RSPG','mean_obj_2RSPG','mean_true_zero_2RSPG','mean_false_zero_2RSPG','best_solution_2RSPG', '-append')
frep=fopen(curr_path2,'a');
fprintf(frep,'     \n');
fprintf(frep,'*** Running 2-RSPG\n');
Results_showing2;
fclose(frep);
count_2RSPG = count_2RSPG + Run_times;
