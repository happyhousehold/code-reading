fprintf('     \n');
if (strcmp(cur_alg, 'RSG') == 1) || (strcmp(cur_alg, '2-RSG') == 1)
    fprintf(frep,'      Iteration limit :%d\n',N_iter);
    
else
    fprintf(frep,'       Number of iterations :%d\n',N_iter);
end
if strcmp(data.problem, 'Least_sqr') == 1
    fprintf(frep,'     \n');
    fprintf(frep,'       Mean of function values :% 5.4f', mean_loss);
    fprintf(frep,'     \n');
    fprintf(frep,'      Variance of function values:% 4.2e\n', var_loss);
elseif strcmp(data.problem, 'Reg') == 1
    fprintf(frep,'     \n');
    fprintf(frep,'       Mean of squared norm of gradient :% 5.4f', mean_grad);
    fprintf(frep,'     \n');
    fprintf(frep,'      Variance of squared norm of gradinet:% 4.2e\n', var_grad);
elseif strcmp(data.problem, 'SVM') == 1
    fprintf(frep,'     \n');
    fprintf(frep,'       Mean of squared norm of gradient :% 5.4f', mean_grad);
    fprintf(frep,'     \n');
    fprintf(frep,'      Variance of squared norm of gradinet:% 4.2e\n', var_grad);
    fprintf(frep,'     \n');
    fprintf(frep,'      Average corresponding missclasification error:% 4.2f\n', mean_loss);
elseif strcmp(data.problem, 'Trig') == 1
    fprintf(frep,'     \n');
    fprintf(frep,'       Mean of squared norm of gradient :% 5.4f', mean_grad);
    fprintf(frep,'     \n');
    fprintf(frep,'      Variance of squared norm of gradinet:% 4.2e\n', var_grad);
    fprintf(frep,'     \n');
    fprintf(frep,'      Mean of corresponding function values:% 4.2e\n', var_grad);
end

