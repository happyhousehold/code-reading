fprintf('  \n');
running = 0;
while running == 0
    mtitle = 'Algorithm Parameters';
    awr = menu(mtitle, ['Number of iterations: ',num2str(N_iter)] , ['Regularization parameter(lambda): ', num2str(lambda)] , ['Evaluation sample size: ', num2str(N_vali)] , ['Number of runs: ', num2str(Run_times)] , ['Magnitude of initial point: ' , num2str(R_ini)],['L_correction: ' , num2str(L_correction)], 'Run Algorithms','Quit');
    if awr ==1
        N_iter = input('Enter number of iterations: ');
    elseif  awr ==2
        lambda = input('Enter the regularization parameter (rho) : ');
        data.lambda = lambda; 
    elseif  awr ==3
        N_vali = input('Enter the validation samples : ');
        data.vali = N_vali;
    elseif  awr ==4
        Run_times = input('Enter run times : ');
    elseif  awr ==5
        R_ini = input('Enter magnitude of initial point: ');
    elseif awr ==6
        L_correction = input('Enter L_correction:  ');
    elseif awr ==7 
        running = 1;
    elseif  awr == 8
        break
    end
end
if  awr == 8
    break
end
