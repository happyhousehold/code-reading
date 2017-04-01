% this file implements a variant of Nesterov's algorithm for solving QP %
function NEST_QP_Original_Spect(control, domain, data, x0)

%%% initialization
start_NEST=clock;
data.L = poweriteration(data.A, domain.m, domain.n, 1e-6);
L = data.L;

%complexity = sqrt(data.L * domain.R^2 / control.epsilon);
%disp(sprintf('complexity: %d', complexity));
n = domain.n;
R = domain.R;
UB = 1e40;
%x0 = reshape(eye(n)*R/n, n^2,1);
xk = x0;
t = 0;
sum_g = 0;
sigma = 0.5;
% Step choose y0 = argmin_{x in X} \frac{L}{sigma} d(x) +
        % \frac{1}{2} [f(x0) + \langle f'(x0), x-x0\rangle]
[f,g]= FirstOrderOracleQP(data,domain,xk);%oracle(data, xk);
UB = f;
yk = ProxmappingSpect(data, domain, sigma/L * g);
sumg = 0;
disp(sprintf('nstep=%d, UB=%.6e \n', 0, full(UB)));
while  t <= control.iter_limit
        % compute zk = argmin sum_g' y + \frac{L}{\sigma} Omega(y): y in X
        sum_g = sum_g + (t+1)/2 * g;
        zk = ProxmappingSpect(data, domain, sigma/L * sum_g);  
        xk = 2/(t+3) * zk + 1/(t+3) * yk;
        [f,g]= FirstOrderOracleQP(data,domain,xk);%oracle(data, xk);
       
        % compute hatx_k = \argmin sigma g(xk) 2 /(L(t+3)) + \left[ Omega(x) - Omega (z_k) -
       % <dOmega(z_k), x - z_k> \right]: x in X
        [V, D] = eig(reshape(zk, n,n)+(1e-10)*eye(n));
       d = diag(1 + log(diag(D)));
       d = reshape(d, n^2,1);        
       hatxk = ProxmappingSpect(data, domain, 2* sigma * g /(L * (t+3)) - d);
      %hatxk = ProjX(data, domain, sigma * g /L - log(zk) -1);
       yk = 2/(t+3) * hatxk + 1 /(t+3) * yk;
       [fyk,gyk]= FirstOrderOracleQP(data,domain,yk);%oracle(data, yk);         
       
       UB = min(UB, min(f, fyk));
        
        if UB < control.epsilon,
            break;
        end;
        
        if mod(t,10) == 0
            disp(sprintf('nstep=%d, UB=%.6e \n', t+1, full(UB)));
        end;
        t = t+1;
end

time_NEST=etime(clock,start_NEST);
disp(sprintf('end of execution.\n'));
str  = sprintf('nstep=%d, UB=%.6e\n,time=%5.2f, L=%5.2f', t, full(UB),time_NEST,data.L);
disp(str);
frep = fopen('report.txt','a');
fprintf(frep,'Summary: %s**\n',str);

fclose(frep);

