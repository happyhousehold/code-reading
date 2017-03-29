%% 2-RSGF-V
if (strcmp(cur_alg, '2-RSGF-V') == 1)
    fprintf('     \n');
    fprintf('*** Running 2-RSGF-V\n');
    num_2RSGF_V = num_2RSGF_V+1;
else
    fprintf('     \n');
    fprintf('*** Running MD-SA-GF\n');
    num_MD_SA_GF = num_MD_SA_GF+1;
end
if N_iter <= 1000
    seed_val = data.seed+120;
else
    seed_val = data.seed+130;
end
Run_number = 0;
total_run = 0;
j=1;
grad_best  = zeros(1, Run_times);
loss_best  = zeros(1, Run_times);
best_solution = zeros(Run_times, data.dim);
num_rnd = S;
while j <= Run_times
    data.vali= ceil(N_iter/(2*S));
    zero_eval_xbar;
    j= j+1;
end
var_grad = var(grad_best);
mean_grad = mean(grad_best);
var_loss = var(loss_best);
mean_loss = mean(loss_best);
zero_Results_showing;
if strcmp(cur_alg, '2-RSGF-V') == 1
    best_solution_2RSGF_V(count_2RSGF_V+1:count_2RSGF_V +Run_times,1) = N_iter;
    best_solution_2RSGF_V(count_2RSGF_V+1:count_2RSGF_V+Run_times,2:data.dim+1) = best_solution;
    var_grad_2RSGF_V(num_2RSGF_V,1) = N_iter;
    var_grad_2RSGF_V(num_2RSGF_V,2) = var_grad;
    mean_grad_2RSGF_V(num_2RSGF_V,1) = N_iter;
    mean_grad_2RSGF_V(num_2RSGF_V,2) = mean_grad;
    var_loss_2RSGF_V(num_2RSGF_V,1) = N_iter;
    var_loss_2RSGF_V(num_2RSGF_V,2) = var_loss;
    mean_loss_2RSGF_V(num_2RSGF_V,1) = N_iter;
    mean_loss_2RSGF_V(num_2RSGF_V,2) = mean_loss;
    save(filename, 'initial_solution','var_grad_2RSGF_V' ,'mean_grad_2RSGF_V','var_loss_2RSGF_V','mean_loss_2RSGF_V','best_solution_2RSGF_V', '-append');
    frep=fopen(curr_path2,'a');
    fprintf(frep,'     \n');
    fprintf(frep,'*** Running 2-RSGF-V\n');
    zero_Results_showing2;
    fclose(frep);
    count_2RSGF_V = count_2RSGF_V + Run_times;
else
    best_solution_MD_SA_GF_GF(count_MD_SA_GF+1:count_MD_SA_GF +Run_times,1) = N_iter;
    best_solution_MD_SA_GF(count_MD_SA_GF+1:count_MD_SA_GF + Run_times,2:data.dim+1) = best_solution;
    var_grad_MD_SA_GF(num_MD_SA_GF,1) = N_iter;
    var_grad_MD_SA_GF(num_MD_SA_GF,2) = var_grad;
    mean_grad_MD_SA_GF(num_MD_SA_GF,1) = N_iter;
    mean_grad_MD_SA_GF(num_MD_SA_GF,2) = mean_grad;
    var_loss_MD_SA_GF(num_MD_SA_GF,1) = N_iter;
    var_loss_MD_SA_GF(num_MD_SA_GF,2) = var_loss;
    mean_loss_MD_SA_GF(num_MD_SA_GF,1) = N_iter;
    mean_loss_MD_SA_GF(num_MD_SA_GF,2) = mean_loss;
    save(filename,'initial_solution','var_grad_MD_SA_GF' ,'mean_grad_MD_SA_GF','var_loss_MD_SA_GF','mean_loss_MD_SA_GF','best_solution_MD_SA_GF', '-append');
    frep=fopen(curr_path2,'a');
    fprintf(frep,'     \n');
    fprintf(frep,'*** Running MD-SA-GF\n');
    zero_Results_showing2;
    fclose(frep);
    count_MD_SA_GF = count_MD_SA_GF + Run_times;
end