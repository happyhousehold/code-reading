% Script for deterministic matrix game
% TODO: Add comments here

%% Init
clear; clc; close all;
addpath(genpath('../../Utilities'));
addpath(genpath('../../Solvers'));

% Parameters
sMOSEK = 'c:\Program Files\mosek\7\toolbox\r2009b'; % MOSEK path
seed = 18; % Seed
m = 1000; % Dimension of y
bSilent = false; % Silent/verbose mode
bSave = true; % Flag for saving all results

addpath(sMOSEK);
RandStream.setDefaultStream(RandStream('mt19937ar','Seed',seed));
if bSave
    sResultDir = sprintf('../../../../Results/D_Game');
    if ~exist(sResultDir, 'dir')
        mkdir(sResultDir);
    end
    sLog = sprintf('%s/log_DG.txt', sResultDir);
    diary(sLog);
    diary on;
end
silent_fprintf(bSilent, '%%-- %s --%% \r\n', datestr(now));
silent_fprintf(bSilent, 'Seed = %g\r\n', seed);

%% Main script
for n = [1000, 10000]
% for n = 1000
    for k = [100, 1000]
%     for k = 1000
        %% Generate problem
        A = randn(k, n);
        K = rand(m, n)*2 - 1;
        Q = A'*A;
        if n>1000 || ~exist('mskqpopt', 'file')
            bPOnly = 1;
        else
            bPOnly = 0;
        end
        LipG = max(max(abs(Q)));
        LipK = max(max(abs(K)));
        silent_fprintf(bSilent, 'm=%g,n=%g,k=%g,LipG=%g,LipK=%g\r\n',...
            m, n, k, LipG, LipK);

        %% General Parameters
        par = [];
        par.bDualityGap = 1;
        par.fhDualityGap = @(x, y)(funG_DualityGap(x, y, K, Q, bPOnly));
        par.bSilent = bSilent;
        par.LipG = LipG;
        par.LipK = LipK;
        par.TolGap = eps;
        
        for MaxIter = [100, 1000, 2000]
%         for MaxIter = 100
            par.MaxIter = MaxIter;
            if k==1000 && MaxIter==2001
                par.OutputInterval = 1;
            else
                par.OutputInterval = MaxIter;
            end
            
            %% Nesterov's smoothing technique
%             silent_fprintf(bSilent, 'Running Nesterov''s algorithm...\r\n');
%             tic;
%             [x, y, etc] = funG_Nesterov(Q, K, par);
%             tEnd = toc;
%             silent_fprintf(bSilent, 'Execution time (sec): %g\r\n', tEnd);
%             xNest = x;
%             yNest = y;
%             etcNest = etc;

            %% AC-SA
%             par.OutputInterval = 1;
            silent_fprintf(bSilent, 'Running AC-SA algorithm...\r\n');
            tic;
            [x, y, etc] = funG_ACSA(Q, K, par);
            tEnd = toc;
            silent_fprintf(bSilent, 'Execution time (sec): %g\r\n', tEnd);
            xACSA = x;
            yACSA = y;
            etcACSA = etc;
            
            %% Nemirovski's prox method
%             silent_fprintf(bSilent, 'Running Nemirovski''s algorithm...\r\n');
%             tic;
%             [x, y, etc] = funG_PM(Q, @(x)(K*x), @(y)((y'*K)'), K, par);
%             tEnd = toc;
%             silent_fprintf(bSilent, 'Execution time (sec): %g\r\n', tEnd);
%             xPM = x;
%             yPM = y;
%             etcPM = etc;

            %% Accelerated primal-dual
%             % The APD function is a universal solver, with more parameters
%             parAPD = par;
%             parAPD.xsize = [n, 1];
%             parAPD.ysize = [m, 1];
%             parAPD.x0 = ones(n, 1)/n;
%             parAPD.y0 = ones(m, 1)/m;
%             parAPD.fhProjx = @(x, dx)(funProxMapEntropy(x, dx));
%             parAPD.fhProjy = @(y, dy)(funProxMapEntropy(y, -dy));
%             parAPD.TolX = eps;
%             parAPD.StepsizePolicy = 1;
%             parAPD.DXYRatio = sqrt(log(n)/log(m));
%             silent_fprintf(bSilent, 'Running APD algorithm...\r\n');
%             tic;
%             [x, y, etc] = funAPD(@(x,t)(Q*x), @(x)(K*x), @(y)((y'*K)'), parAPD);
%             tEnd = toc;
%             silent_fprintf(bSilent, 'Execution time (sec): %g\r\n', tEnd);
%             xAPD = x;
%             yAPD = y;
%             etcAPD = etc;

            %% Linearized version of Chambolle-Pock
%             parPD = parAPD;
%             parPD.StepsizePolicy = 0;
%             silent_fprintf(bSilent, 'Running linearized version of Chambolle-Pock algorithm...\r\n');
%             tic;
%             [x, y, etc] = funAPD(@(x,t)(Q*x), @(x)(K*x), @(y)((y'*K)'), parPD);
%             tEnd = toc;
%             silent_fprintf(bSilent, 'Execution time (sec): %g\r\n', tEnd);
%             xPD = x;
%             yPD = y;
%             etcPD = etc;

            %% Plot duality gap
            if k==1000 && MaxIter==2001
                h = figure;
                if bPOnly
                    plot(repmat(1:MaxIter, 4, 1)', ...
                        [etcNest.PrimalObjectiveValue, etcPM.PrimalObjectiveValue, etcPD.PrimalObjectiveValue, etcAPD.PrimalObjectiveValue]);
                    legend('Nesterov', 'Nemirovski', 'Linearized Chambolle-Pock', 'APD');
                    xlabel('Iteration');
                    ylabel('Primal objective value');
                    ylim([0, .01]);
                else
                    plot(repmat(1:MaxIter, 4, 1)', ...
                        [etcNest.DualityGap, etcPM.DualityGap, etcPD.DualityGap, etcAPD.DualityGap]);
                    legend('Nesterov', 'Nemirovski', 'Linearized Chambolle-Pock', 'APD');
                    xlabel('Iteration');
                    ylabel('Duality gap');
                    ylim([0, .01]);
                end
                if bSave
                    print(h, '-dpdf', ...
                        sprintf('%s/DG_%g_m%gn%gk%gN%g', ...
                        sResultDir, datenum(date), m, n, k, MaxIter));
                end
            end
        end
        fprintf('\r\n');
    end
end

diary off;