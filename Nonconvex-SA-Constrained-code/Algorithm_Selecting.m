mtitle = 'Algorithm Selecting';
method_selection = menu (mtitle, 'RSPG', '2-RSPG', '2-RSPG-V','Quit');
if method_selection == 1
    cur_alg = 'RSPG';
    cur_alg2 = 'RSPG';
    Input_Parameters;
elseif  method_selection == 2
    cur_alg = '2-RSPG';
    cur_alg2 = 'RSPG_2';
    Input_Parameters;
elseif  method_selection == 3
    cur_alg = '2-RSPG-V';
    cur_alg2 = 'RSPG_2V';
    Input_Parameters;
elseif  method_selection == 4
    indicator = 1;
    break;
end
