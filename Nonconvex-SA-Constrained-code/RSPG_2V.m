%% Hybrid RSPG
m = data.m;
%N_iter = ceil(N_total/m);
if (strcmp(cur_alg, '2-RSPG-V') == 1)
    fprintf('     \n');
    fprintf('*** Running 2-RSPG-V\n');
    num_2RSPG_V = num_2RSPG_V + 1;
end
if N_total <= 1000
    seed_val = data.seed+260;
elseif   N_total <= 5000
    seed_val = data.seed+280;
else
    seed_val = data.seed+300;
end
Run_number = 0;
j=1;
grad1_best  = zeros(1, Run_times);
grad2_best  = zeros(1, Run_times);
loss_best  = zeros(1, Run_times);
obj_best  = zeros(1, Run_times);
true_zero=zeros(1, Run_times);
false_zero=zeros(1, Run_times);
best_solution = zeros(Run_times, data.dim);
num_rnd = S;
while j <= Run_times
    data.vali= ceil(N_iter/(2*S));
    eval_xbar;
    j= j+1;
end
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
best_solution_2RSPG_V(count_2RSPG_V+1:count_2RSPG_V + Run_times,1) = N_iter;
best_solution_2RSPG_V(count_2RSPG_V+1:count_2RSPG_V + Run_times,2:data.dim+1) = best_solution;
var_grad1_2RSPG_V(num_2RSPG_V,1) = N_iter;
var_grad1_2RSPG_V(num_2RSPG_V,2) = var_grad1;
var_grad2_2RSPG_V(num_2RSPG_V,1) = N_iter;
var_grad2_2RSPG_V(num_2RSPG_V,2) = var_grad2;
mean_grad1_2RSPG_V(num_2RSPG_V,1) = N_iter;
mean_grad1_2RSPG_V(num_2RSPG_V,2) = mean_grad1;
mean_grad2_2RSPG_V(num_2RSPG_V,1) = N_iter;
mean_grad2_2RSPG_V(num_2RSPG_V,2) = mean_grad2;
var_loss_2RSPG_V(num_2RSPG_V,1) = N_iter;
var_loss_2RSPG_V(num_2RSPG_V,2) = var_loss;
mean_loss_2RSPG_V(num_2RSPG_V,1) = N_iter;
mean_loss_2RSPG_V(num_2RSPG_V,2) = mean_loss;
var_obj_2RSPG_V(num_2RSPG_V,1) = N_iter;
var_obj_2RSPG_V(num_2RSPG_V,2) = var_obj;
mean_obj_2RSPG_V(num_2RSPG_V,1) = N_iter;
mean_obj_2RSPG_V(num_2RSPG_V,2) = mean_obj;
mean_true_zero_2RSPG_V(num_2RSPG_V,1) = N_iter;
mean_true_zero_2RSPG_V(num_2RSPG_V,2) = mean_true_zero;
mean_false_zero_2RSPG_V(num_2RSPG_V,1) = N_iter;
mean_false_zero_2RSPG_V(num_2RSPG_V,2) = mean_false_zero;
save(filename,'var_grad1_2RSPG_V' ,'mean_grad1_2RSPG_V','var_grad2_2RSPG_V' ,'mean_grad2_2RSPG_V','var_loss_2RSPG_V','mean_loss_2RSPG_V','var_obj_2RSPG_V','mean_obj_2RSPG_V','mean_true_zero_2RSPG_V','mean_false_zero_2RSPG_V','best_solution_2RSPG_V', '-append');
frep=fopen(curr_path2,'a');
fprintf(frep,'     \n');
fprintf(frep,'*** Running 2-RSPG-V\n');
Results_showing2;
fclose(frep);
count_2RSPG_V = count_2RSPG_V + Run_times;

