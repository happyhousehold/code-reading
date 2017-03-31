% Main Program
clc
clear
problem_select;
if  strcmp(data.problem, 'S3VM') == 1
    initializer;
    General_sample_generation;
    fprintf(frep,'********** Instance: %s ',filename);
    fprintf(frep,'     \n');
    fclose(frep);
    filename = [curr_path,'\Results-',filename];
    save(filename,'data');
    indicator = 0;
    Estimating_Parameters;
    common_prog;
end
