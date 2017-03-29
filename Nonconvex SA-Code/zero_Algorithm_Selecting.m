mtitle = 'Algorithm Selecting';
method_selection = menu (mtitle, 'RSGF', '2-RSGF', '2-RSGF-V', 'MD-SA-GF','Change the problem' ,'Quit');
if method_selection == 1
    cur_alg = 'RSGF';
    cur_alg2 = 'RSGF';
    zero_Input_Parameters;
elseif  method_selection == 2
    cur_alg = '2-RSGF';
    cur_alg2 = 'RSGF_2';
    zero_Input_Parameters;
elseif  method_selection == 3
    cur_alg = '2-RSGF-V';
    cur_alg2 = 'RSGF_2V';
    zero_Input_Parameters;
elseif  method_selection == 4
    cur_alg = 'MD-SA-GF';
    cur_alg2 = 'RSGF_2V';
    zero_Input_Parameters;
elseif  method_selection == 5
    Main_Program;
elseif  method_selection == 6
    indicator = 1;
    break;
end
