% this file implements Nesterov's smoothing algorithm for solving 
% the eigenvalue problem
function NEST_Eigen(control, domain, data)

%%% initialization
%complexity = sqrt(data.L * domain.R^2 / control.epsilon);
%disp(sprintf('complexity: %d', complexity));

% first estimate the operator norm 
% A_norm = norm(data.A0, 'fro');
% for i = 1:data.n,
%     L1 = norm(data.A0+data.A{i}, 'fro');
%     if A_norm < L1,
%         A_norm = L1;
%     end
% end

A_norm = normest(abs(data.A0));
for i = 1:data.n,
    L1 = normest(abs(data.A0+data.A{i}));
    if A_norm < L1,
        A_norm = L1;
    end
end

% second compute the smoothing parameter
%mu = control.epsilon / (2 * log(domain.m));

D1 = log(domain.n);
D2 = log(domain.m);
sigma1 = 1;
sigma2 = 1;
mu = 2 * A_norm / (control.iter_limit+1) * sqrt(D1 / (sigma1 * sigma2 * D2));

data.mu = mu; %this value will be used in the oracle

% estimate the Lipschitz constant
L = A_norm^2 / mu;

L = L * control.stepfactor;

% choose x0 as the minimizer of the entropy function
xk = ones(domain.n,1)/domain.n;
omega = -log(domain.n);
domega = zeros(domain.n,1);
sigma = 1.0;

% choose y0 = argmin { L w(x) / sigma + g(x0)' x / 2}
%           = argmin { w(x) + sigma g(x0)' x / (2L) }
[f,g,f0] = Oracle(data, xk, 0);
sum_g = g / 2.0;
xi = sigma * g / (2 * L);
[yk,omeganew,domeganew]=ProjMapping(xi,domain);

UB = 1e9;
UB0 = 1e9;
k = 0;
    
while  k <= control.iter_limit
    % Find zk = argmin { L w(x) / sigma + sum_g * x: x \in Q}
    % zk = argmin{w(x) + sigma * sum_g * x / L: x \in Q}
    xi = sigma * sum_g / L;
    [zk,omeganew,domeganew]=ProjMapping(xi,domain);

    tauk = 2.0 / (k+3);
    xk = tauk * zk  + (1-tauk)* yk;        

    [f,g,f0]= oracle(data, xk, 0);
    sum_g = sum_g + (k+2) * g / 2.0;

    %hxk = V_Q(zk, sigma tauk g(xk) / L)
    %    = argmin <sigma tauk g(xk) / L, x> + w(x) - <w'(zk), x>
    xi = sigma * tauk * g / L - domeganew;
    [hxk,omeganew,domeganew]=ProjMapping(xi,domain);

    yk = tauk * hxk + (1 - tauk) * yk;

    [fyk,gyk,f0yk]= oracle(data, xk, 0);

    if f < UB,
        x_ubt = xk;
        UB = f;
        UB0 = f0;
    end
    
    if fyk < UB,
        x_ubt = yk;
        UB = fyk;
        UB0=f0yk;
    end

 %   if mod(k,100) == 0
        disp(sprintf('nstep=%d, UB=%.6e UB0=%.6e\n', k, UB, UB0));
 %   end;
    k = k+1;
end
    
disp(sprintf('end of execution.\n'));
str  = sprintf('nstep=%d, UB=%.6e, UB0=%.6e', k, UB, UB0);
disp(str);
frep = fopen('report.txt','a');
fprintf(frep,'Summary: %s**\n',str);
fclose(frep);

