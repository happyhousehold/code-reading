%% Choosing the problem
mtitle = 'Choose the problem you want to solve:';
problem_selection = menu (mtitle, 'Least Square', 'Non-convex Regression', 'Non-convex supprot vector machine', 'Trigonometric function','Inventory','Quit');
if problem_selection == 1
    data.problem = 'Least_sqr';
elseif  problem_selection == 2
    data.problem ='Reg';
elseif  problem_selection == 3
    data.problem ='SVM';
elseif problem_selection == 4
    data.problem ='Trig';
elseif  problem_selection == 5
    data.problem ='Invt';
elseif  problem_selection == 6
    data.problem ='';
end
