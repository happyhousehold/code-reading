if strcmp(cur_alg, 'RSGF') == 1
    N_description = 'Iteration limit: ';
else
    N_description = 'Number of iterations: ';
end
awr = menu(mtitle, [N_description ,num2str(N_iter)] ,['initial_solution: [', num2str(initial_solution'),']'], ['Evaluation sample size (K): ', num2str(N_vali)] , ['Number of runs: ', num2str(Run_times)] , ['Run:', cur_alg], 'Change the algorithm');
if awr == 1
    N_iter = input('Enter the iteration limit: ');
elseif  awr == 2
    ini_sol = input('Enter the initial solution : ','s');
    initial_solution = str2num(ini_sol)';
elseif  awr == 3
    N_vali = input('Enter the size of evaluation sample (K): ');
    data.vali = N_vali;
elseif  awr == 4
    Run_times = input('Enter number of runs : ');
elseif  awr == 5
    eval(cur_alg2);
elseif  awr == 6
    running = 1;
    break
end
