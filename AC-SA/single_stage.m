fprintf('     \n');
fprintf('*** Running Single-stage AC-SA\n');
for j = 1: Run_times
    tic
    z = initial_solution;
    z_ag = z;
    t = 1;
    while  t <= N_iter
        % Step1
        alfa = 2/(t+1);
        beta = 4 * L / (t * (t+1));
        cx_p1 = (1-alfa)*(mu + beta) /( beta + (1 - alfa^2) * mu);
        z_md = cx_p1 * z_ag + (1.0-cx_p1) *z;
        cx_p2 = alfa * mu / (mu + beta);
        z_ps = cx_p2 * z_md + (1.0 - cx_p2) * z;
        index = t;
        Data_generation;
        SGradient;
        z = (-alfa/(mu + beta))*GR + z_ps;
        z_ag = alfa*z + (1-alfa)*z_ag;
        t = t+1;
    end
    final_solution = z_ag';
    objective_computation;
end
Computing_Data;
% file_name = [filename,'.txt'];
% frep=fopen(file_name,'a');
frep=fopen('Results.txt','a');
count_alg = count_alg +1;
if count_alg == 1
    fprintf(frep,'     \n');
    fprintf(frep,'**********Data: %s , rho: %6.5g\n',filename, lamda);
end
fprintf(frep,'     \n');
fprintf(frep,'*** Running Single-stage AC-SA\n');
fprintf(frep,'     \n');
fprintf(frep,'       Runtime:% 4.2f     Final objective value:% 6.2f\n', final_data);
fclose(frep);
