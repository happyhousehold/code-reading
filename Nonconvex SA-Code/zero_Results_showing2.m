fprintf('     \n');
if (strcmp(cur_alg, 'RSGF') == 1) || (strcmp(cur_alg, '2-RSGF') == 1)
    fprintf(frep,'      Iteration limit :%d\n',N_iter);
    
else
    fprintf(frep,'       Number of iterations :%d\n',N_iter);
end
if strcmp(data.problem, 'Invt') == 1
    fprintf(frep,'     \n');
    fprintf(frep,'       Mean of squared norm of gradient :% 5.4f', mean_grad);
    fprintf(frep,'     \n');
    fprintf(frep,'      Variance of squared norm of gradinet:% 4.2e\n', var_grad);
    fprintf(frep,'     \n');
    fprintf(frep,'       Mean of inventory cost :% 5.4f', mean_loss);
    fprintf(frep,'     \n');
end
