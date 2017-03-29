%% Choosing the problem
mtitle = 'Choose the problem you want to solve:';
problem_selection = menu (mtitle, 'Convex Least Square problem', 'Non-convex supprot vector machine problem ', 'Simulation-based inventory optimization problem','Quit');
if problem_selection == 1
    data.problem = 'Least_sqr';
elseif  problem_selection == 2
    data.problem ='SVM';
elseif problem_selection == 3
    data.problem ='Invt';
elseif  problem_selection == 4
    data.problem ='';
end
