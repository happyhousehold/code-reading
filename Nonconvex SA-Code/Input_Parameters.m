%% Setting algorithm parameters
running = 0;
while running == 0
    mtitle = 'Algorithm Parameters';
    if  (strcmp(data.problem, 'Reg') == 1) || (strcmp(data.problem, 'SVM') == 1)
        if (strcmp(cur_alg, '2-RSG') == 1) || (strcmp(cur_alg, '2-RSG-V') == 1)
            Input1;
        else
            Input2;
        end
    else
        if (strcmp(cur_alg, '2-RSG') == 1) || (strcmp(cur_alg, '2-RSG-V') == 1)
            Input3;
        else
            Input4;
        end
    end
end
common_prog;
