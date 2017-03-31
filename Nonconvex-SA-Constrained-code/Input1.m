if strcmp(cur_alg, '2-RSPG') == 1
    N_description = 'Iteration limit: ';
else
    N_description = 'Number of iterations: ';
end
awr = menu(mtitle, [N_description,num2str(N_iter)],['Mini-batch size: ',num2str(data.m)] ,['# of Candidate solutions: ',num2str(S)] ,['Regularization parameters [lambda1,lambda2,lambda3]: [', num2str(data.lambda1),',',num2str(data.lambda2),',',num2str(data.lambda3),']'] , ['Evaluation sample size (K): ', num2str(N_vali)] , ['Number of runs: ', num2str(Run_times)] , ['Run: ', cur_alg],'Change the algorithm');
if awr == 1
    N_iter = input('Enter the iteration limit: ');
elseif  awr == 2
    data.m = input('Enter the mini-batch size: ');
elseif  awr == 3
    S = input('Enter the number of candidate solutions: ');
elseif  awr == 4
    lambda_vector = input('Enter the regularization parameters [lambda1,lambda2,lambda3] : ');
    data.lambda1= lambda_vector(1);
    data.lambda2= lambda_vector(2);
    data.lambda3= lambda_vector(3);
    Estimating_Parameters;
elseif  awr == 5
    N_vali = input('Enter the size of evaluation sample (K): ');
    data.vali = N_vali;
elseif  awr == 6
    Run_times = input('Enter number of runs : ');
elseif  awr == 7
    eval(cur_alg2);
elseif  awr == 8
    running = 1;
    break
end
