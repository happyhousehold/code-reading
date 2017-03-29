% Main Program
clc
clear
problem_selecting;
curr_path = pwd;
curr_path = [curr_path,'\Result'];
curr_path2 = [curr_path,'\Results.txt'];
frep = fopen(curr_path2,'a');
fprintf(frep,'     \n');
if strcmp(data.problem, '') ~= 1
    if (strcmp(data.problem, 'Invt') ~= 1)
        initializer;
        General_sample_generation;
        fprintf(frep,'********** Instance: %s ',filename);
        fprintf(frep,'     \n');
        fclose(frep);
        indicator = 0;
        filename = [curr_path,'\Results-',filename];
        save(filename,'data');
        Estimating_Parameters;
        common_prog;
    else
        zero_initializer;
        fprintf(frep,'********** Inventory problem');
        fprintf(frep,'     \n');
        fclose(frep);
        indicator = 0;
        filename = [curr_path,'\Results-Inventory'];
        save(filename,'data');
        zero_Estimating_Parameters;
        zero_common_prog;
    end
end