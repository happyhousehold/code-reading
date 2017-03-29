z = initial_solution;
if strcmp(cur_alg2, 'RSG_2V') == 1
    r = N_iter;
else
    generating_rnd_num;
    r = r1;
end
Run_number = Run_number + 1;
classic_SA;
final_sol(Run_number,:) = z;
final_solution = final_sol(Run_number,:)';

