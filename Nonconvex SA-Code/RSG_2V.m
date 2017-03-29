%% Hybrid RSG
if (strcmp(cur_alg, '2-RSG-V') == 1)
    fprintf('     \n');
    fprintf('*** Running 2-RSG-V\n');
    num_2RSG_V = num_2RSG_V + 1;
else
    fprintf('     \n');
    fprintf('*** Running MD-SA\n');
    num_MD_SA = num_MD_SA+1;
end
if N_iter <= 1000
    seed_val = data.seed+260;
elseif   N_iter <= 5000
    seed_val = data.seed+280;
else
    seed_val = data.seed+300;
end
Run_number = 0;
j=1;
grad_best  = zeros(1, Run_times);
loss_best  = zeros(1, Run_times);
best_solution = zeros(Run_times, data.dim);
num_rnd = S;
while j <= Run_times
    data.vali= N_iter/(2*S);
    eval_xbar;
    j= j+1;
end
var_grad = var(grad_best);
mean_grad = mean(grad_best);
var_loss = var(loss_best);
mean_loss = mean(loss_best);
Results_showing;
if strcmp(cur_alg, '2-RSG-V') == 1
    best_solution_2RSG_V(count_2RSG_V+1:count_2RSG_V + Run_times,1) = N_iter;
    best_solution_2RSG_V(count_2RSG_V+1:count_2RSG_V + Run_times,2:data.dim+1) = best_solution;
    var_grad_2RSG_V(num_2RSG_V,1) = N_iter;
    var_grad_2RSG_V(num_2RSG_V,2) = var_grad;
    mean_grad_2RSG_V(num_2RSG_V,1) = N_iter;
    mean_grad_2RSG_V(num_2RSG_V,2) = mean_grad;
    var_loss_2RSG_V(num_2RSG_V,1) = N_iter;
    var_loss_2RSG_V(num_2RSG_V,2) = var_loss;
    mean_loss_2RSG_V(num_2RSG_V,1) = N_iter;
    mean_loss_2RSG_V(num_2RSG_V,2) = mean_loss;
    save(filename,'var_grad_2RSG_V' ,'mean_grad_2RSG_V','var_loss_2RSG_V','mean_loss_2RSG_V','best_solution_2RSG_V', '-append');
    frep=fopen(curr_path2,'a');
    fprintf(frep,'     \n');
    fprintf(frep,'*** Running 2-RSG-V\n');
    Results_showing2;
    fclose(frep);
    count_2RSG_V = count_2RSG_V + Run_times;
else
    best_solution_MD_SA(count_MD_SA+1:count_MD_SA +Run_times,1) = N_iter;
    best_solution_MD_SA(count_MD_SA+1:count_MD_SA+Run_times,2:data.dim+1) = best_solution;
    var_grad_MD_SA(num_MD_SA,1) = N_iter;
    var_grad_MD_SA(num_MD_SA,2) = var_grad;
    mean_grad_MD_SA(num_MD_SA,1) = N_iter;
    mean_grad_MD_SA(num_MD_SA,2) = mean_grad;
    var_loss_MD_SA(num_MD_SA,1) = N_iter;
    var_loss_MD_SA(num_MD_SA,2) = var_loss;
    mean_loss_MD_SA(num_MD_SA,1) = N_iter;
    mean_loss_MD_SA(num_MD_SA,2) = mean_loss;
    save(filename,'var_grad_MD_SA' ,'mean_grad_MD_SA','var_loss_MD_SA','mean_loss_MD_SA','best_solution_MD_SA', '-append');
    frep=fopen(curr_path2,'a');
    fprintf(frep,'     \n');
    fprintf(frep,'*** Running MD-SA\n');
    Results_showing2;
    fclose(frep);
    count_MD_SA = count_MD_SA + Run_times;
end
