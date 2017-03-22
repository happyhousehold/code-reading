% basic_run.m

clear;clc;
root=pwd;root = root(1:end-8);
addpath([root '/problems'],[root '/solvers'],[root '/solvers/utilities']);
    
% set up the run
n = 2^10;
delta = 0.3;            % m/n, m = round(delta*n)
rho = 0.2;              % k/m, k = round(rho*m)
Ameth = 0;              % see getData.m for codes
xmeth = 0;              %  "      "      "    "
opts = fpc_opts([]);    % see fpc_opts.m for options
sig1 = 1E-2;            % std. dev. of signal noise
sig2 = 1E-8;            %  "    "   "  measurement noise
full = true;            % whether to use a full or approximate M matrix
mu = [];                % mu to use--[] means recommended
sig1est = sig1;         % estimate of sig1 used by getM_mu
sig2est = sig2;         %    "     "  sig2  "   "    "
alpha = 0.5;            % parameter for chi^2 value
nseMult = 3;            % in debias: nse = nseMult*sigma
seed = 0;               % if nonempty, will be seed for rand and randn
opts.fullMu = false;    % if false, just update mu as opts.eta*mu
opts.scale = false;     % if false, fpc assumes problem already scaled
printPlots = false;     % make pdfs of plots?

% plots
signals = false;
xsConvergence = true;
paretoPlot = false;
paretoPlotPlusProgress = false;

% problem size
m = round(delta*n);
k = round(rho*m);

% get problem
data_t = cputime;
[A,b,xs,xsn] = getData(m,n,k,Ameth,xmeth,sig1,sig2,seed);
data_t = cputime - data_t;
disp([num2str(data_t),' s to get the problem.']);

Mmu_t = cputime;
[M,mu_rec,A,b,sig,kap,tau,M12] = getM_mu(full,m,n,Ameth,A,b,sig1est,sig2est,alpha);
if isempty(mu), mu = mu_rec; end
if ~isempty(M), A = M12*A; b = M12*b; M = []; end
if ~isempty(tau), opts.tau = tau; end
opts.kappa = kap;
Mmu_t = cputime - Mmu_t;
disp([num2str(Mmu_t),' s to estimate M and mu.']);

opts.xs = xs;

% fpc-basic
solve_t = cputime;
Out = fpc(n,A,b,mu,M,opts);
if Ameth == 5, Out.x = real(Out.x); end
solve_t = cputime - solve_t;
disp([num2str(Out.itr),' iterations and ',num2str(solve_t),...
    ' s to solve the problem to rel. err. ',num2str(Out.n2re(end),'%5.3g'),...
    ' with basic fpc.']);

% de-bias the solution
db_t = cputime;
nse = nseMult*sig;
x = debias(m,n,Out.x,A,b,M,nse);
db_t = cputime - db_t;
n2re_db = norm(x - xs)/norm(xs);
disp([num2str(db_t),' s to debias.  Resulting rel. err. is ',...
    num2str(n2re_db),'.']);

% plots
if signals
    figure('PaperSize',[4.5,3.75],'PaperPosition',[0,0,4.5,3.75]);
    subplot('Position',[0.12,0.78,0.85,0.20]);
    plot([1:n]',xs,'k-','LineWidth',0.75);
    set(gca,'XLim',[0,n],'YLim',[-4,4],'XTickLabel',[],'YTick',[-4,0,4]);
    ylabel('x_s');
    subplot('Position',[0.12,0.56,0.85,0.20]);
    plot([1:n]',xsn,'k-','LineWidth',0.75);
    set(gca,'XLim',[0,n],'YLim',[-4,4],'XTickLabel',[],'YTick',[-4,0,4]);
    ylabel('x_s + \epsilon_1');
    subplot('Position',[0.12,0.34,0.85,0.20]);
    plot([1:n]',Out.x,'k-','LineWidth',0.75);
    set(gca,'XLim',[0,n],'YLim',[-4,4],'XTickLabel',[],'YTick',[-4,0,4]);
    ylabel('x');
    subplot('Position',[0.12,0.12,0.85,0.20]);
    plot([1:n]',x,'k-','LineWidth',0.75);
    set(gca,'XLim',[0,n],'YLim',[-4,4],'YTick',[-4,0,4]);
    ylabel('Debiased');
    xlabel('Component');
    if printPlots, print nse_rec -dpdf; end
end

if xsConvergence
    figure('PaperSize',[4,4],'PaperPosition',[0,0,4,4]);
    semilogy([0:Out.itr],Out.n2re);
    xlabel('Iteration');
    ylabel('||x - x_s||/||x_s||');
    if printPlots, print errVitr -dpdf; end
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
        figure('PaperSize',[4,4],'PaperPosition',[0,0,4,4]);
        plot(lamS,phiS,'k-'); xlabel('||x||_1'); ylabel('||Ax - b||_M'); 
        if printPlots, print pareto -dpdf; end
    end
    if paretoPlotPlusProgress
        figure('PaperSize',[4,4],'PaperPosition',[0,0,4,4]);
        plot(lamS,phiS,'k-'); xlabel('||x||_1'); ylabel('||Ax - b||_M'); 
        hold on; plot(Out.lam,phi,'b-'); legend('pareto curve','fpc progress');
        if printPlots, print paretoWprogress -dpdf; end
    end
end