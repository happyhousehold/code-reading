min_grad = 1e+10;
best_i = 0;
if strcmp(cur_alg2, 'RSGF_2V') == 1
    generating_rnd_num;
    gh = r1;
    smooth_zero_SA;
    if strcmp(cur_alg, '2-RSGF-V') == 1
        for i = 1:S
            final_solution = all_iterates (gh(i),:)';
            zero_evaluation;
            if gradient_final < min_grad
                min_grad = gradient_final;
                best_i = gh(i);
            end
        end
        final_solution = all_iterates (best_i,:)';
    else
        average_sol = mean(all_iterates);
        final_solution = average_sol';
    end
    seed_val = seed_val+1;
else
    for i = 1: S
        smooth_zero_SA;
        zero_evaluation;
        if gradient_final < min_grad
            min_grad = gradient_final;
            best_i = Run_number;
        end
        seed_val = seed_val+1;
    end
    final_solution = final_sol(best_i,:)';
end
Run_number = 0;
clear final_sol
clear all_iterates
data.vali = N_vali;
zero_evaluation;
best_solution(j,:) = final_solution';
grad_best(j) = gradient_final;
loss_best(j) = final_loss;
