
if (strcmp(data.problem, 'Least_sqr') == 1) || (strcmp(data.problem, 'Reg') == 1)
    %% Regression problem
    er = normrnd(0,data.st);
    x = sprand(1, data.dim, data.spr);
    y = dot(x, z_sample) + er;
    
elseif (strcmp(data.problem, 'SVM') == 1)
    %% Non-convex SVM problem
    x = ceil(sprand(data.dim, 1, data.spr));
    y = sign(dot(x, z_sample));
    
end
