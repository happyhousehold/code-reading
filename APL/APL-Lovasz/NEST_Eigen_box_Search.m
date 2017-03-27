% this file implements Nesterov's algorithm for solving Compsite problem 
% with line search
function NEST_Eigen_box_Search(control, domain, data)

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
data.L = A_norm^2 / mu;

UB = 1e9;
x0 = zeros(domain.n,1);
xk = x0;
vk = xk;
t = 0;
sum_g = zeros(domain.n,1);

A = 0;
L = 1;
while  t <= control.iter_limit
        % Step1
        
      %% line search for Lk
      while 1 == 1
        alpha = 1 / L + sqrt( 1 / L^2 + 2 * A / L);
        yk = (A * xk + alpha * vk) / (A + alpha);
        
        [f,g,phi]= oracle(data, yk);
      
%        Tyk = CompositeProj(g - L*yk, L, domain.n);
        
        Tyk = yk - g / L;
        Tyk = max(0,Tyk);
        Tyk = min(Tyk,data.DomainType);
        
        [fTy, gTy, phiTy] = oracle(data, Tyk);
        
        if L >= data.L,
            break;
        end;
        
        sTy = sign(Tyk);
        g_gTy = gTy + sTy;
      
        if g_gTy' * (yk - Tyk) < g_gTy'*g_gTy / L,
            L = L *2;
        else
            break;
        end;
      end
      
      M = L;
      L = M/2;
      xk = Tyk;
      [fxk,gxk,phixk]= oracle(data, xk);
      sum_g = sum_g + alpha * gxk;
      A = A + alpha;
  
      % compute vk = argmin sum_g' y / A +  \|y - x0\|^2 / (2 A) + \|y\|_1
      vk = x0 - sum_g;
      vk = max(0,vk);
      vk = min(0,data.DomainType);
        
      UB = min(UB, phixk);
        
        if mod(t,100) == 0
            disp(sprintf('nstep=%d, UB=%.6e\n', t, UB));
        end;
        t = t+1;
end

disp(sprintf('end of execution.\n'));
str  = sprintf('nstep=%d, UB=%.6e, L=%.6e\n', t, UB, L);
disp(str);
frep = fopen('report.txt','a');
fprintf(frep,'Summary: %s**\n',str);

fclose(frep);

