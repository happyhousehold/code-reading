% this file implements a variant of Nesterov's algorithm for solving QP %
function NEST_QP_Original_Spect_Eu(control, domain, data, x0)

%%% initialization
start_NEST=clock;
data.L = poweriteration(data.A, domain.m, domain.n, 1e-2)
%opts.issym = 1;
%opts.disp = 0;
%[V,D] = eigs(sparse(data.A'*data.A), 1, 'la', opts);
%data.L = D;
    
%complexity = sqrt(data.L * domain.R^2 / control.epsilon);
%disp(sprintf('complexity: %d', complexity));
data.L = data.L;
n=domain.n;
R=domain.R;
UB = 1e9;
%x0 = reshape(eye(n)*R/n, n^2,1);
xk = x0;
t = 0;
sum_g = 0;
type = domain.type;

    while  t <= control.iter_limit
        % Step1
        
        [f,g]= FirstOrderOracleQP(data,domain,xk);%oracle(data, xk);
        % compute y_k = argmin g' y + L \| y - x_k\|^2/2: \|y\|_\infty \le R
        %yk = xk - g / data.L;        
        %yk = max(min(yk, domain.R),0);
        
        yk = ProjSpectahedron(data,domain,g/data.L - xk);
                
        [fyk,gyk]= FirstOrderOracleQP(data,domain,yk);%oracle(data, yk);  
        
        sum_g = sum_g + (t+1) * g / 2.0;
        
        % compute zk = argmin sum_g' y + L \|y\|^2 / 2: \|y\|_\infty \le R
        %zk = -sum_g / data.L;
        %zk = max(min(zk, domain.R),0);
        zk = ProjSpectahedron(data,domain,sum_g/data.L);        
           
        xk = 2 * zk / (t+3) + (t+1) * yk / (t+3);
        %z_md = (1-alfa)*z_ag + alfa*z;
        
        UB = min(UB, min(f, fyk));
        
        if UB < control.epsilon,
            break;
        end;
        
        if mod(t,10) == 0
            disp(sprintf('nstep=%d, UB=%.6e \n', t, full(UB)));
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

