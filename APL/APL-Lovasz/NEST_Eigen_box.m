% this file implements Nesterov's algorithm for solving Compsite problem %
function NEST_Eigen_box(control, domain, data)

%%% initialization

% first estimate the operator norm 
A_norm = norm(data.A0, 'fro');
for i = 1:data.n,
    L1 = norm(data.A0+data.A{i}, 'fro');
    if A_norm < L1,
        A_norm = L1;
    end
end

% second compute the smoothing parameter
%mu = control.epsilon / (2 * log(domain.m));
D1 = domain.n;
D2 = log(domain.m);
sigma1 = 1;
sigma2 = 1;
mu = 2 * A_norm / (control.iter_limit+1) * sqrt(D1 / (sigma1 * sigma2 * D2));

data.mu = mu; %this value will be used in the oracle

% estimate the Lipschitz constant
L = A_norm^2 / mu;

UB = 1e9;
xk = zeros(domain.n,1);
t = 0;
A = 0;
sum_g = 0;

    while  t <= control.iter_limit
        % Step1
        
        [f,g,phi]= oracle(data, xk);
        % compute y_k = argmin g' y + L \| y - x_k\|^2/2
        % yk = CompositeProj(g - data.L*xk, data.L, domain.n);
        yk = xk - g / L;
        yk = max(0,yk);
        yk = min(yk,data.DomainType);
        
        [fyk,gyk,phiyk]= oracle(data, yk);  
        
        sum_g = sum_g + (t+1) * g / 2.0;
        A = A + (t+1) / 2.0;
        
        % compute zk = argmin sum_g' y / A + L \|y\|^2 / (2A)
        %zk = CompositeProj(sum_g / A, data.L / A, domain.n);
        zk = - sum_g / L;
        zk = max(0, zk);
        zk = min(zk,data.DomainType);
        
        xk = 2 * zk / (t+3) + (t+1) * yk / (t+3);
        %z_md = (1-alfa)*z_ag + alfa*z;
        
        UB = min(UB, min(phi, phiyk));
        
%        if UB-data.opt < control.epsilon,
%            break;
%        end;
        
 %       if mod(t,100) == 0
            disp(sprintf('nstep=%d, UB=%.6e \n', t, UB));
 %       end;
        t = t+1;
    end
    
disp(sprintf('end of execution.\n'));
str  = sprintf('nstep=%d, UB=%.6e, dataL=%.6e\n', t, UB, L);
disp(str);
frep = fopen('report.txt','a');
fprintf(frep,'Summary: %s**\n',str);

fclose(frep);

