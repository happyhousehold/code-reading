% Demo for random sampling matrix game

addpath(genpath('../../Utilities'));
addpath(genpath('../../Solvers'));
addpath('../LS_TV'); % For ClA_Operator
RandStream.setDefaultStream(RandStream('mt19937ar','Seed',18));

% m = 100;
% n = 100;
m = 10000;
n = 10000;
k = 100;
nRun = 100;
bSave = true;

if bSave
    sResultDir = sprintf('../../../../Results/S_Game');
    if ~exist(sResultDir, 'dir')
        mkdir(sResultDir);
    end
    sLog = sprintf('%s/log_SG.txt', sResultDir);
    diary(sLog);
    diary on;
end

[cind, rind] = meshgrid(1:m, 1:n);

for iType = [1, 2]
% for iType = 1
    for alpha = [2, 1, .5]
%     for alpha = 2        
        %% Generate matrices K and Q
        disp('Generating Q...');
        A = randn(k, n);
        At = A';
        Q = At*A;
        LipG = max(max(abs(Q)));
        objQ = ClA_operator(@(x)(At*(A*x)), []);
        disp('Generating K...');
        switch iType
            case 1
                % % K of Type I
                K = ((rind + cind - 1)/(2*n - 1)).^alpha;
            case 2
                % K of Type II
                K = ((abs(rind-cind) + 1)/(2*n-1)).^alpha;
        end
        LipK = max(max(abs(K)));
        
        % Stochastic gradients
        Kt = K';
        fhK = @(x)(K(:, randsample(n,1,true,x)));
        fhKt = @(y)(Kt(:, randsample(m,1,true,y)));
        
        fhDualityGap = @(x, y)(funG_DualityGap(x, y, K, objQ, 1));
        
        state = get(RandStream.getDefaultStream, 'state');
        
        %% Initial POBJ
        x0 = ones(n, 1)./n; y0 = ones(m, 1)./m;
        [~, POBJ0] = fhDualityGap(x0, y0);
        fprintf('Type %d,alpha=%g,LipG=%g,LipK=%g,Inital POBJ: %g\r\n', iType, alpha, LipG, LipK, POBJ0);
        
        %% General parameters
        par = [];
        par.bSilent = true;
        par.bDualityGap = 1;
        par.fhDualityGap = fhDualityGap;
        par.LipK = LipK;
        par.LipG = LipG;
        
        for MaxIter = [100, 1000, 2000]
%         for MaxIter = 1000
            fprintf('N=%d\n', MaxIter);
            par.MaxIter = MaxIter;
            par.OutputInterval = MaxIter;
            
            %% Robust SA
            set(RandStream.getDefaultStream, 'state', state);
            parSA = par;
            parSA.M_star = sqrt(2*log(n)*(LipG + LipK)^2 + 2*log(m)*LipK^2);
            POBJ_SA = zeros(nRun, 1);
            disp('Running robust SA...');
            tic;
            for i = 1:nRun
                [~, ~, etc] = funG_MDSA(objQ, fhK, fhKt, K, parSA);
                POBJ_SA(i) = etc.PrimalObjectiveValue(end);
            end
            fprintf('Robust SA, Mean: %g, Std: %g\n', mean(POBJ_SA), std(POBJ_SA));
            tEnd = toc;
            fprintf('Avg. Execution time (sec): %g\n', tEnd/nRun);
            
            %% SMP
            set(RandStream.getDefaultStream, 'state', state);
            parSMP = par;
            parSMP.M = LipK * sqrt(2 * (log(n) + log(m)));
            POBJ_SMP = zeros(nRun, 1);
            disp('Running SMP...');
            tic;
            for i = 1:nRun
                [~, ~, etc] = funG_PM(objQ, fhK, fhKt, K, par);
                POBJ_SMP(i) = etc.PrimalObjectiveValue(end);
            end
            fprintf('SMP, Nean: %g, Std: %g\n', mean(POBJ_SMP), std(POBJ_SMP));
            tEnd = toc;
            fprintf('Avg. Execution time (sec): %g\n', tEnd/nRun);
            
            %% Accelerated primal-dual
            set(RandStream.getDefaultStream, 'state', state);
            POBJ_APD = zeros(nRun, 1);
            disp('Running APD...');
            tic;
            % The APD function is a universal solver, with more parameters
            parAPD = par;
            parAPD.xsize = [n, 1];
            parAPD.ysize = [m, 1];
            parAPD.x0 = ones(n, 1)/n;
            parAPD.y0 = ones(m, 1)/m;
            parAPD.fhProjx = @(x, dx)(funProxMapEntropy(x, dx));
            parAPD.fhProjy = @(y, dy)(funProxMapEntropy(y, -dy));
            parAPD.TolX = eps;
            parAPD.StepsizePolicy = 1;
            parAPD.DXYRatio = sqrt(log(n)/log(m));
            parAPD.sigma_xDX = LipK / sqrt(2*log(n));
            parAPD.sigma_yDY = LipK / sqrt(2*log(m));
            for i = 1:nRun
                [~, ~, etc] = funAPD(@(x,t)(At*(A*x)), fhK, fhKt, parAPD);
                POBJ_APD(i) = etc.PrimalObjectiveValue(end);
            end
            fprintf('S-APD, Mean: %g, Std: %g\n', mean(POBJ_APD), std(POBJ_APD));
            tEnd = toc;
            fprintf('Avg. Execution time (sec): %g\n', tEnd/nRun);
            
        end
    end
end

diary off;