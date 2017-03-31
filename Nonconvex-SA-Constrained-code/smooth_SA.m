z = initial_solution;
%z = (abs(initial_solution)>=threshold_val).*initial_solution;
%z_0 = (z==0)+z;
if (strcmp(cur_alg2, 'RSPG_2V') == 1)
    r = N_iter;
else
    generating_rnd_num;
    r = r1;
end
Run_number = Run_number + 1;
classic_SA;
final_sol(Run_number,:) = z;
final_solution = final_sol(Run_number,:)';



