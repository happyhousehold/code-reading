mtitle = 'Algorithm Selecting';
method_selection = menu (mtitle, 'Single-stage AC-SA', 'Multi-stage AC-SA', 'Shrinking single-stage AC-SA' , 'Shrinking multi-stage AC-SA','Classic SA' , 'Batch-learning', 'Change algorithm parameters','Change initial sapmle' ,'Quit');
if method_selection == 1
    single_stage;
elseif  method_selection == 2
    multi_stage;
elseif method_selection == 3
    single_stage_shrinkage;
elseif  method_selection == 4
    multi_stage_shrinkage;
elseif method_selection == 5
    classic_SA;
elseif  method_selection == 6
    deterministic;
elseif  method_selection == 7
    fprintf('  \n');
    Input_Parameters;
    Estimating_Parameters;
elseif  method_selection == 8
    running = 0;
    Sample_generation;
    Estimating_Parameters;
elseif  method_selection == 9
    indicator = 1;
end