min_grad = 1e+10;
best1_i =0;
min_loss = 1e+10;
best2_i =0;
if strcmp(cur_alg, 'HBRSGF') == 1
    gholam = randi(N_iter,1, S);
    smooth_zero_SA;
    for i = 1:S
        final_solution = all_iterates (gholam(i),:)';
        zero_evaluation;
        if gradient_final < min_grad
            min_grad = gradient_final;
            best1_i = gholam(i);
        end
        if final_loss < min_loss
            min_loss = final_loss;
            best2_i = gholam(i);
        end
    end
    final_solution = all_iterates (best1_i,:)';
    Run_number=0;
    data.vali = N_vali;
    zero_evaluation;
    grad_best1(j) = gradient_final;
    loss_best1(j) = final_loss;
    final_solution = all_iterates (best2_i,:)';
    data.vali = N_vali;
    zero_evaluation;
    grad_best2(j) = gradient_final;
    loss_best2(j) = final_loss;
    clear final_sol
    clear all_iterates
end


