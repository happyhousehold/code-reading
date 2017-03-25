fprintf('     \n');
fprintf('*** Running Classic SA\n');
for j = 1 : Run_times
    tic
    % c is the strong convexity parameter
    c = 10*min(eig(A));
    
    teta = 1/c;
    z = initial_solution;
    t = 1;
    while  t <= N_iter
        index = t;
        Data_generation;
        z_md = z;
        SGradient;
        gamma = teta/t;
        z = z - gamma*GR;
        t = t+1;
    end
    final_solution = z';
    objective_computation;
end
Computing_Data;
% file_name = [filename,'.txt'];
% frep=fopen(file_name,'a');
frep=fopen('Results.txt','a');
count_alg = count_alg +1;
if count_alg == 1
    fprintf(frep,'     \n');
    fprintf(frep,'********** Data: %s , rho: %6.5g\n',filename, lamda);
end
fprintf(frep,'     \n');
fprintf(frep,'*** Running Classic SA\n');
fprintf(frep,'     \n');
fprintf(frep,'       Runtime:% 4.2f     Final objective value:% 6.2f\n', final_data);
fclose(frep);

