% one_run.m

clear;
    
% set up the run
n = 2^12;
delta = 0.1;            % m/n, m = round(delta*n)
rho = 0.1;              % k/m, k = round(rho*m)
Ameth = 6;              % see getData.m for codes
xmeth = 0;              %  "      "      "    "
opts = fpc_opts([]);    % see fpc_opts.m for options
sig1 = 1E-3;            % std. dev. of signal noise
sig2 = 1E-3;            %  "    "   "  measurement noise
full = false;           % whether to use a full or approximate M matrix
opts.fullMu = false;    % if false, just update mu as opts.eta*mu
mu = [];                % mu to use--[] means recommended
sig1est = sig1;         % estimate of sig1 used by getM_mu
sig2est = sig2;         %    "     "  sig2  "   "    "
alpha = 0.5;            % parameter for chi^2 value
nseMult = 3;            % in debias: nse = nseMult*sigma

% plots
xsConvergence = true;
paretoPlot = true;
paretoPlotPlusProgress = true;

% problem size
m = round(delta*n);
k = round(rho*m);

% get problem
data_t = cputime;
[A,b,xs,xsn,picks] = getData(m,n,k,Ameth,xmeth,sig1,sig2,1978);
data_t = cputime - data_t;
disp([num2str(data_t),' s to get the problem.']);

Mmu_t = cputime;
[M,mu,A,b,sig,kap,tau,M12] = getM_mu(full,mu,m,n,Ameth,A,b,sig1est,sig2est,alpha);
if ~isempty(M), A = M12*A; b = M12*b; M = []; end
if ~isempty(tau), opts.tau = tau; end
opts.kappa = kap;
Mmu_t = cputime - Mmu_t;
disp([num2str(Mmu_t),' s to estimate M and mu.']);

opts.xs = xs;

% fpc-basic
solve_t = cputime;
Out = fpc(n,A,b,mu,M,opts,picks);
solve_t = cputime - solve_t;
disp([num2str(Out.itr),' iterations and ',num2str(solve_t),...
    ' s to solve the problem to rel. err. ',num2str(Out.n2re(end),'%5.3g'),...
    ' w/ fpc-basic.']);

% de-bias the solution
db_t = cputime;
nse = nseMult*sig;
x = debias(m,n,Out.x,A,b,M,nse,picks);
db_t = cputime - db_t;
n2re_db = norm(x - xs)/norm(xs);
disp([num2str(db_t),' s to de-bias.  Resulting rel. err. is ',...
    num2str(n2re_db),'.']);

% plots
if xsConvergence
    figure;
    semilogy([0:Out.itr],Out.n2re);
    xlabel('iteration');
    ylabel('||x - x_s||/||x_s||');
end

if paretoPlot || paretoPlotPlusProgress
    phi = []; phiS = []; lamS = [];
    start = 1;
    for i = 1:length(Out.itrs)
        last = sum(Out.itrs(1:i)) + 1;
        phi = [phi; sqrt((2/Out.mus(i))*(Out.f(start:last) - Out.lam(start:last)))];
        if i < length(Out.itrs)
            phi(end) = phi(end)*sqrt(Out.mus(i)/Out.mus(i+1));
        end
        phiS = [phiS; phi(end)]; lamS = [lamS; Out.lam(last)];
        start = last + 1;
    end
    if paretoPlot
        figure;
        plot(lamS,phiS,'k-'); xlabel('||x||_1'); ylabel('||Ax - b||_M'); 
    end
    if paretoPlotPlusProgress
        figure;
        plot(lamS,phiS,'k-'); xlabel('||x||_1'); ylabel('||Ax - b||_M'); 
        hold on; plot(Out.lam,phi,'b-');
    end
end