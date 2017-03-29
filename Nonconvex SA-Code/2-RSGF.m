%% 2-RSGF
cur_alg = '2-RSGF';
j = 1;
grad_best  = zeros(1, Run_times);
loss_best  = zeros(1, Run_times);
num_rnd = 1;
while j <= Run_times
    data.vali = N_iter/2;
    zero_eval_xbar;
    j = j+1;
end
var_grad_2RSGF(Run_count) = var(grad_best);
mean_grad_2RSGF(Run_count) = mean(grad_best);
var_loss_2RSGF (Run_count) = var(loss_best);
mean_loss_2RSGF(Run_count) = mean(loss_best);
grad_best  = zeros(1, Run_times);
loss_best  = zeros(1, Run_times);