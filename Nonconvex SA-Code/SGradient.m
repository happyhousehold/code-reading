
if (strcmp(data.problem, 'Least_sqr') == 1) || (strcmp(data.problem, 'Reg') == 1)
    %% Non-convex Regression
    Data_generation;
    GR = 2*(dot(z, x)-y)*x'+ data.lambda*p*z.*(z.^2+ep^2*ones(d,1)).^(p/2-1);

elseif (strcmp(data.problem, 'SVM') == 1)
    %% Non-convex SVM PROBLEM
    Data_generation;
    GR = y*(tanh(y*dot(z,x))^2-1)*x+ 2*data.lambda*z;
    
elseif (strcmp(data.problem, 'Trig') == 1)
    %% Trigonometric functions
    GR_0 = (data.dim - sum(cos(z)))*ones(data.dim, 1)+cumsum(ones(data.dim, 1)).*(ones(data.dim, 1)-cos(z))-sin(z);
    GR = 2*(sum(GR_0)*sin(z)+GR_0.*(cumsum(ones(data.dim, 1)).*sin(z)-cos(z)))+normrnd(0, data.st, data.dim,1);
    
end