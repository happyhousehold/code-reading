min_grad = 1e+6;
min_loss  =1e+6;
best_i =0;
if strcmp(cur_alg2, 'RSG_2V') == 1
    generating_rnd_num;
    gh = r1;
    smooth_SA;
    if strcmp(cur_alg, '2-RSG-V') == 1
        for i = 1:S
            final_solution = all_iterates (gh(i),:)';
            objective_computation;
            Computing_Data;
            if strcmp(data.problem, 'Least_sqr') == 1
                if final_loss < min_loss
                    min_loss = final_loss;
                    best_i = gh(i);
                end
            else
                if gradient_final < min_grad
                    min_grad = gradient_final;
                    best_i = gh(i);
                end
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
        smooth_SA;
        objective_computation;
        Computing_Data;
        if strcmp(data.problem, 'Least_sqr') == 1
            if final_loss < min_loss
                min_loss = final_loss;
                best_i = Run_number;
            end
        else
            if gradient_final < min_grad
                min_grad = gradient_final;
                best_i = Run_number;
            end
        end
        seed_val = seed_val+1;
    end
    final_solution = final_sol(best_i,:)';
end
Run_number = 0;
clear final_sol
clear all_iterates
data.vali = N_vali;
objective_computation;
Computing_Data;
best_solution(j,:) = final_solution';
grad_best(j) = gradient_final;
loss_best(j) = final_loss;
