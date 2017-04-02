% TODO: Add comments

%% Init
clear; 
clc; 
close all;
addpath(genpath('../../Solvers'));
addpath(genpath('../../Utilities'));

% Parameters
sigma = .001; % Standard deviation of noise
seed = 18; % Seed
bSilent = false; % Silent/verbose mode
bSave = 0; % Flag for saving all results

RandStream.setDefaultStream(RandStream('mt19937ar','Seed',seed));
if bSave
    sResultDir = sprintf('../../../../Results/LS_TV');
    if ~exist(sResultDir, 'dir')
        mkdir(sResultDir);
    end
    sLog = sprintf('%s/log_LS_TV.txt', sResultDir);
    diary(sLog);
    diary on;
end
silent_fprintf(bSilent, '%%-- %s --%% \r\n', datestr(now));
silent_fprintf(bSilent, 'Seed = %g\r\n', seed);

% for sInstance = {'PPI', 'Gaussian'}
for sInstance = {'PPI'}
    %% Generate data
    switch sInstance{:}
        case 'PPI'
            load('../../../../Data/PPI/data1.mat');
            xTrue = u0;
            [nRow, nCol, nCh] = size(sense_map);
            opA = @(x)(bsxfun(@times, fft2(bsxfun(@times, x, sense_map)), p));
            Noise = sigma*sqrt(nRow*nCol/2)*(randn(nRow, nCol, nCh) + 1i*randn(nRow, nCol, nCh));
            b = bsxfun(@times, opA(xTrue) + Noise, p);
            opAt = @(y)(sum(ifft2(bsxfun(@times, y, p)) .* conj(sense_map), 3));
            objA = ClA_operator(opA, opAt);
            LipG = max(max(abs(sum(sense_map, 3))))^2;
%             lwReg = [1e-3, 1e-4]/(nRow*nCol);
%             lwReg = 10.^(-13:2);
            lwReg = [2e-5:1e-5:5e-5,9e-6:-1e-6:5e-6]; 
%             lwReg = 1e-5;
        case 'Gaussian'
            xTrue = phantom(64);
            [nRow, nCol] = size(xTrue);
            nSample = ceil(nRow*nCol/2);
            A = randn(nSample, nRow*nCol) / sqrt(nRow*nCol);
            b = A * xTrue(:) + randn(nSample, 1) * sigma;
            LipG = eigs(A' * A, 1);
            objA = ClA_operator(@(x)(A*x(:)), @(b)(reshape(A'*b, [nRow, nCol])));
            lwReg = [1e-2, 1e-3];
    end
    %% Necessary function handles, constants and parameters for the solver
    % Gradient of the quadratic term
    fhGradG = @(x) (objA' * (objA * x - b));

    % Parameters for the uniform solver
    par = [];
    par.xsize = [nRow, nCol];
    par.wsize = [nRow, nCol, 2];
    par.ysize = [nRow, nCol, 2];
    par.LipG = LipG;
    par.xTrue = xTrue;
    par.TolX = eps;
    par.bSilent = bSilent;
    par.bObjectiveValue = true;
    par.bRelativeError = true;
    par.fhRelativeError = @funRelativeL2Error;
    par.fhProjy = @funProxMapEuclL21;
    
    for wReg = lwReg
        par.LipK = sqrt(8) * wReg;
        % Operators K, Kt
        fhK = @(x)(funTVGrad(x, wReg, 0));
        fhKt = @(y)(funTVNegDiv(y, wReg, 0));
        % Function handle for calculating energies
        fhPOBJ = @(x)(L2TVEnergy(objA, x, b, fhK));
        par.fhObjectiveValue = fhPOBJ;

        silent_fprintf(bSilent, 'Instance %s,lambda=%g\r\nLipG=%g,LipK=%g\r\n',...
            sInstance{:}, wReg, par.LipG, par.LipK);

%         for MaxIter = [50 100 150 1000]
        for MaxIter=1000
            par.MaxIter = MaxIter;
%             par.OutputInterval = MaxIter;
            par.OutputInterval = 150;
            par.bPlot = 1;
            
            %% Accelerated Primal Dual (APD)
            par.StepsizePolicy = 1;
            silent_fprintf(bSilent, 'Running APD...\r\n');
            tic;
            [x, y, etc] = funAPD(fhGradG, fhK, fhKt, par);
            tEnd = toc;
            silent_fprintf(bSilent, 'Execution time (sec): %g\r\n', tEnd);
            xAPD = x;
            yAPD = y;
            etcAPD = etc;
%             %% APD, unbounded version (APD-U)
%             par.StepsizePolicy = 2;
%             silent_fprintf(bSilent, 'Running APD, unbounded version...\r\n');
%             tic;
%             [x, y, etc] = funAPD(fhGradG, fhK, fhKt, par);
%             tEnd = toc;
%             silent_fprintf(bSilent, 'Execution time (sec): %g\r\n', tEnd);
%             xAPDU = x;
%             yAPDU = y;
%             etcAPDU = etc;
%             %% Linearized version of Chambolle-Pock (LPD)
%             par.StepsizePolicy = 0;
%             par.bErgodic = 0;
%             silent_fprintf(bSilent, 'Running linearized version of Chambolle-Pock algorithm (LPD)...\r\n');
%             tic;
%             [x, y, etc] = funAPD(fhGradG, fhK, fhKt, par);
%             tEnd = toc;
%             silent_fprintf(bSilent, 'Execution time (sec): %g\r\n', tEnd);
%             xLPD = x;
%             yLPD = y;
%             etcLPD = etc;
%             %% Linearized version of Chambolle-Pock (LPD), ergodic solution
%             par.StepsizePolicy = 0;
%             par.bErgodic = 1;
%             silent_fprintf(bSilent, 'Running linearized version of Chambolle-Pock algorithm (LPD), ergodic solution...\r\n');
%             tic;
%             [x, y, etc] = funAPD(fhGradG, fhK, fhKt, par);
%             tEnd = toc;
%             silent_fprintf(bSilent, 'Execution time (sec): %g\r\n', tEnd);
%             xLPD = x;
%             yLPD = y;
%             etcLPD = etc;
        end
    end
    silent_fprintf(bSilent, '\r\n');
end

diary off;