function [f,g,f0,diff_eig] = Oracle_simple(data, x)
% note that f_mu = mu ln [sum_i exp (lambda_i(x) /mu)] - mu log m
% Assemble the matrix
%ExA = data.A0;
%for i=1:data.n,
%    ExA = ExA + data.A{i} * x(i);
%end

longvector = data.exA * [1;x];
sumA = reshape(longvector, data.m, data.m);
% Compute the eigenvalue decomposition
options.disp = 0;
[V,D] = eigs(sumA,1,'LA', options);
eig_val1 = D;

% Compute the difference between the first and second largest eigenvalue
diff_eig = 0;

f0 = eig_val1;
g = zeros(data.n,1);

    V1 = V; % the eigevector corresponding to the largest eigenvalue    
    gouter0 = V1 * V1';
    gouter = reshape(gouter0, [],1);
    g = data.exA(:,2:data.n+1)' * gouter;
    f = f0;



