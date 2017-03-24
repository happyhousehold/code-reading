% fprintf('  \n');
% running = 0;
% while running == 0
% mtitle = 'Algorithm Parameters';
% awr = menu(mtitle, ['Number of iterations: ',num2str(N_iter)] , ['Regularization parameter(rho): ', num2str(lamda)] , ['Validation samples: ', num2str(N_vali)] , ['Run times: ', num2str(Run_times)] , ['Magnitude initial point: ' , num2str(R_ini)], 'Run Algorithms','Quit');
%     if awr ==1
%         N_iter = input('Enter number of iterations: ');
%     elseif  awr ==2
%         lamda = input('Enter the regularization parameter (rho) : ');
%     elseif  awr ==3
%         N_vali = input('Enter the validation samples : ');
%     elseif  awr ==4
%         Run_times = input('Enter run times : ');
%     elseif  awr ==5
%         R_ini = input('Enter magnitude initial point: ');
%     elseif awr ==6
%         running = 1;
%     elseif  awr == 7
%         return;
%     end
% end


% N_iter = input('Enter number of iterations: ');
% lamda = input('Enter the regularization parameter (rho) : ');
% N_vali = input('Enter the validation samples : ');
% Run_times = input('Enter run times : ');
% R_ini = input('Enter magnitude initial point: ');


mtitle = 'Algorithm Parameters';
awr = menu(mtitle, ['Number of iterations: ',num2str(N_iter)] , ['Regularization parameter(rho): ', num2str(lamda)] , ['Validation samples: ', num2str(N_vali)] , ['Run times: ', num2str(Run_times)] , ['Magnitude initial point: ' , num2str(R_ini)], 'Run Algorithms','Quit');
