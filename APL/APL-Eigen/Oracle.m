function [f,g,f0,diff_eig,lb] = Oracle(data, x, mode)
%%% if mode = 0, then we report both gradient and function values
%%% if mode = 1, then we only report the function for the original
%%% problem
% note that f_mu = mu ln [sum_i exp (lambda_i(x) /mu)] - mu log m
% Assemble the matrix
%ExA = data.A0;
%for i=1:data.n,
%    ExA = ExA + data.A{i} * x(i);
%end

f = -inf;
g = zeros(data.n,1);
f0 = -inf;
diff_eig = 0;
lb = -inf;

longvector = data.exA * [1;x];
sumA = reshape(longvector, data.m, data.m);

% if mode = 1, we only need to compute the original funciton value,
% so only need to compute the max eigenvalue
options.disp = 0;
if mode == 1   
    [V,D] = eigs(sumA,1,'LA', options);
    f0 = D;
else
    if data.mu == 0 % no smoothing
        [V,D] = eigs(sumA,1,'LA', options);       
        f0 = D;
        V1 = V; % the eigevector corresponding to the largest eigenvalue    
        gouter0 = V1 * V1';
        gouter = reshape(gouter0, [],1);
        g = data.exA(:,2:data.n+1)' * gouter;
        f = f0;
    else
        [V,D] = eig(full(sumA));
        eig_val = diag(D);
        % Compute the difference between the first and second largest eigenvalue
        % diff_eig = eig_val(data.m) - eig_val(data.m-1);
        mx = max(eig_val);
        f0 = mx;
        
        sm_mx = max(eig_val / data.mu);
        normalized = exp(eig_val/data.mu -sm_mx);
        sum_normalized = sum(normalized);
        new_eig = normalized / sum_normalized;

        gouter0 = V * diag(new_eig) * V';
        gouter = reshape(gouter0, [],1);
        g = data.exA(:,2:data.n+1)' * gouter;
        f = data.mu * (sm_mx + log(sum_normalized));
        f = f - data.mu * log(data.m);
    end
end

%%% to compute a lower bound on the optimal value
%%% note that a dual solution is given by 
if mode == 0    
    lb = data.exA(:,1)' * gouter;
    if data.DomainType == 0, % standard simplex
        lb = lb+min(gouter);
    else %box
        lb = lb + min(0,gouter);
    end
end
% % % Compute the eigenvalue decomposition
% % [V,D] = eig(full(sumA));
% % eig_val = diag(D);
% % % Compute the difference between the first and second largest eigenvalue
% % diff_eig = eig_val(data.m) - eig_val(data.m-1);
% % mx = max(eig_val);
% % f0 = mx;
% % f = -inf;
% % g = zeros(data.n,1);
% % 
% % if mode == 0, %%% compute gradient   
% %     if data.mu ==0 %no smoothing
% %         V1 = V(:,data.m); % the eigevector corresponding to the largest eigenvalue    
% %         gouter0 = V1 * V1';
% %         gouter = reshape(gouter0, [],1);
% %         g = data.exA(:,2:data.n+1)' * gouter;
% %         f = f0;
% %     else    % with smoothing
% %         sm_mx = max(eig_val / data.mu);
% %         normalized = exp(eig_val/data.mu -sm_mx);
% %         sum_normalized = sum(normalized);
% %         new_eig = normalized / sum_normalized;
% % 
% %         gouter0 = V * diag(new_eig) * V';
% %         gouter = reshape(gouter0, [],1);
% %         g = data.exA(:,2:data.n+1)' * gouter;
% %         %g = zeros(data.n, 1);
% %         %for i = 1:data.n,
% %         %    g(i) = trace(gouter * data.A{i});
% %         %end
% %         f = data.mu * (sm_mx + log(sum_normalized));
% %         f = f - data.mu * log(data.m);
% %     end
% % end
% % 
