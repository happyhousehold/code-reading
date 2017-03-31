fprintf('     \n');
if (strcmp(cur_alg, 'RSPG') == 1) || (strcmp(cur_alg, '2-RSPG') == 1)
    fprintf(frep,'      Iteration limit :%d\n',N_iter);
    
else
    fprintf(frep,'       Number of iterations :%d\n',N_iter);
end
fprintf('      Mini-batch size :%d\n',data.m);
fprintf(frep,'     \n');
fprintf(frep,'       Mean of squared norm of gradient :% 5.4f', mean_grad1);
fprintf(frep,'     \n');
fprintf(frep,'      Variance of squared norm of gradinet:% 4.2e\n', var_grad1);
fprintf(frep,'     \n');
fprintf(frep,'      Average corresponding objective values:% 6.4f\n', mean_obj);


