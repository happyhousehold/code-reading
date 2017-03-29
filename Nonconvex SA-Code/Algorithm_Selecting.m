mtitle = 'Algorithm Selecting';
method_selection = menu (mtitle, 'RSG', '2-RSG', '2-RSG-V', 'MD-SA','Change the problem' ,'Quit');
if method_selection == 1
    cur_alg = 'RSG';
    cur_alg2 = 'RSG';
    Input_Parameters;
elseif  method_selection == 2
    cur_alg = '2-RSG';
    cur_alg2 = 'RSG_2';
    Input_Parameters;
elseif  method_selection == 3
    cur_alg = '2-RSG-V';
    cur_alg2 = 'RSG_2V';
    Input_Parameters;
elseif  method_selection == 4
    cur_alg = 'MD-SA-GF';
    cur_alg2 = 'RSG_2V';
    Input_Parameters;
elseif  method_selection == 5
    Main_Program;
elseif  method_selection == 6
    indicator = 1;
    break;
end
