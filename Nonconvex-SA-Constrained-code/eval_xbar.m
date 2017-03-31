min_grad = 1e+6;
min_loss  =1e+6;
best_i =0;
if (strcmp(cur_alg2, 'RSPG_2V') == 1)
    generating_rnd_num;
    gh = r1;
    smooth_SA;
    if (strcmp(cur_alg, '2-RSPG-V') == 1)
        for i = 1:S
            final_solution = all_iterates (gh(i),:)';
            objective_computation;
            Computing_Data;
            if gradient_smooth < min_grad
                min_grad = gradient_smooth;
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
        smooth_SA;
        objective_computation;
        Computing_Data;
        if gradient_smooth < min_grad
            min_grad = gradient_smooth;
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
objective_computation2;
Computing_Data;
best_solution(j,:) = final_solution';
grad1_best(j) = gradient_original;
grad2_best(j) = gradient_smooth;
%loss_best(j) = final_loss;
obj_best(j) = objective_value;
final_solution_trc = (abs(final_solution)>=data.threshold).*final_solution;
qq=find(data.sample==0);
true_zero(j)=sum(final_solution_trc(qq)==0);
false_zero(j)=sum(final_solution_trc==0)-true_zero(j);





