fprintf('     \n');
fprintf('*** Running Shrinking multi-stage AC-SA\n');
variance_estimate = var(x_initial')+ mean(x_initial').^2;
sigma = 2*st*sqrt(max(variance_estimate));
for j = 1 : Run_times
    tic
    p = initial_solution;
    lb_mu = 2;
    ek = (1/N)*(y_initial - x_initial*p)'*(y_initial - x_initial*p) + lamda*norm(p)^2;
    conf_level = 0.01;
    K_big = ceil(log(max(1 , 2*ek/epsln)));
    lamda_K_big = -log(conf_level/(2*K_big));
    k = 1;
    N_sum = 0;
    index = 0 ;
    while N_sum < N_iter
        Nk = ceil(max(sqrt(L/mu) , 5*(M^2+(sigma^2)*max(1 , 2*lamda_K_big^2))/(mu*ek)));
        gamma = max(2*L, sqrt(((M^2+sigma^2)*mu/(3*ek))*(Nk+1)*(Nk+2)*Nk));
        Rk = 5*sqrt (4*ek/lb_mu);
        if Nk > (N_iter - N_sum)
            Nk = N_iter - N_sum;
        end
        AC_SA_sub;
        ek = ek/2;
        N_sum = N_sum + Nk;
        k = k+1;
    end
    final_solution = p';
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
fprintf(frep,'*** Running Shrinking multi-stage AC-SA\n');
fprintf(frep,'     \n');
fprintf(frep,'       Runtime:% 4.2f     Final objective value:% 6.2f\n', final_data);
fclose(frep);


