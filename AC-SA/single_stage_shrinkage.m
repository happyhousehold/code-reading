fprintf('     \n');
fprintf('*** Running Shrinking single-stage AC-SA\n');
variance_estimate = var(x_initial')+ mean(x_initial').^2;
sigma = 2*st*sqrt(max(variance_estimate));
for j = 1 : Run_times
    tic
    %%% setup parameters for the shrinking AC-SA algorithm
    p = initial_solution;
    ek = (1/N)*(y_initial - x_initial*p)'*(y_initial - x_initial*p) + lamda*norm(p)^2;
    lb_mu = 2;
    Rk = sqrt (2*ek/mu);
    conf_level = 0.1;
    lambda_t = log ( N_iter / conf_level);
    S = max(4 * sqrt(2 * L / mu), 64 * lambda_t * sigma^2 / (mu * ek) );
    S = max(S, (sigma^2 / (mu * ek)) * (32 *sqrt(14) * lambda_t / 3)^2);
    z = initial_solution;
    z_ag = z;
    t = 1;
    sh_ind = 1;
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
        beta_p = mu + beta;
        z_est = (-alfa/beta_p)*GR + z_ps;
        if norm(z_est - p) <= Rk
            z = z_est;
        else
            c = norm(beta_p*(z_ps - p)-alfa*GR);
            z = (Rk/c)*(beta_p*(z_ps - p)- alfa*GR) + p;
        end
        z_ag = alfa*z + (1-alfa)*z_ag;
        %%% shrinking the feasible set
        tmp_ind = ceil(log(t / S));
        if tmp_ind >= sh_ind
            p = z_ps;
            Rk = Rk / sqrt(2);
            sh_ind = sh_ind + 1;
        end
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
    fprintf(frep,'********** Data: %s , rho: %6.5g\n',filename, lamda);
end
fprintf(frep,'     \n');
fprintf(frep,'*** Running Shrinking single-stage AC-SA\n');
fprintf(frep,'     \n');
fprintf(frep,'       Runtime:% 4.2f     Final objective value:% 6.2f\n', final_data);
fclose(frep);

