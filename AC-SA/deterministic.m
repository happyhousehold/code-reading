fprintf('     \n');
fprintf('*** Running Batch-learning\n');
for j = 1 : Run_times
    tic
    er = normrnd(0,st,1,N_iter)';
    x = rand(N_iter,d);
    y = x*z_sample + er;
    B = inv(N*lamda*eye(d)+x'*x);
    final_solution = (B*x'*y)';
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
fprintf(frep,'*** Running Batch-learning\n');
fprintf(frep,'     \n');
fprintf(frep,'       Runtime:% 4.2f     Final objective value:% 6.2f\n', final_data);
fclose(frep);
