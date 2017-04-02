% Demo for model sum_{i=1}^{nCH}|MFS_jx-f_j|^2/2 + wReg*TV(x), where M is
% the mask, F is the 2D Fourier transform matrix, S_j's are sensitivity
% maps, and f_j's are k-space datasets. The solution x of the model is the
% reconstructed image from partially parallel MRI.

%% Init
clear; 
clc; close all;

% General parameter
seed = 18;
sigma = 1e-3;
wReg = 1e1;
MaxIter = 300;
LipGRelaxation = 1;
LipKRelaxation = 1;
% lAlg = {'APD', 'LPD', 'BOS', 'USL', 'BOSVS'}; % Algorithms to run
% lAlg = {'APD', 'USL'}; 
lAlg = { 'APD'}; 
%% Load data
addpath(genpath('../../Solvers'));
addpath(genpath('../../Utilities'));
addpath(genpath('../../../External/SBB_and_BOSVS'));
addpath(genpath('../../../External/mftv'));
RandStream.setDefaultStream(RandStream('mt19937ar','Seed',seed));

% Generate PPI observation
% load('../../../Data/PPI/data1.mat');
load data2.mat
% [m, n, k] = size(sense_map);
% opA = @(x)(bsxfun(@times, fft2(bsxfun(@times, x, sense_map)), p)); % Matrix A
% opAt = @(y)(sum(ifft2(bsxfun(@times, y, p)) .* conj(sense_map), 3)); % Matrix A'/(m*n)
% 
% A = ClA_operator(opA, opAt);
% knoise = sigma*(randn(m,n,k) + 1i*randn(m,n,k))/sqrt(2);
% f = bsxfun(@times, opA(u0) + knoise, p);


xTrue= u0(100:199,100:199);
[nRow, nCol] = size(xTrue);
nSample = ceil(nRow*nCol/2);
A = randn(nSample, nRow*nCol);
b = A * xTrue(:) + randn(nSample, 1) * sigma;
% LipG = eigs(A' * A, 1);
LipG =2.90e4;



%% Necessary function handles, constants and parameters for the uniform solver
% Gradient of the fidelity term
% fhGradG = @(x) (A' * (A * x - f));
% % Operator K, Kt
% fhK = @(x)(funTVGrad(x, wReg, 0));
% fhKt = @(y)(funTVNegDiv(y, wReg, 0));
% % Function handle for calculating energies
% fhEnergy = @(x)(.5*sum(sum(sum(abs(A*x-f).^2)))/(m*n) + sum(sum(sqrt(sum(abs(fhK(x)).^2,3)))));
% % Energy of ground truth
% TrueEnergy = fhEnergy(u0);

  objA = ClA_operator(@(x)(A*x(:)), @(b)(reshape(A'*b, [nRow, nCol])));
  fhGradG = @(x) (objA' * (objA * x - b));
  fhK = @(x)(funTVGrad(x, wReg, 0));
  fhKt = @(y)(funTVNegDiv(y, wReg, 0));
  fhEnergy = @(x)(.5*norm(A*x(:)-b)^2 + sum(sum(sqrt(sum(abs(fhK(x)).^2,3)))));
  TrueEnergy = fhEnergy(xTrue);
  
  
% Parameters for the uniform solver
par = [];
LipK = sqrt(8) * wReg;
    par.xsize = [nRow, nCol];
    par.wsize = [nRow, nCol, 2];
    par.ysize = [nRow, nCol, 2];
par.LipG = LipG * LipGRelaxation;
par.LipK = LipK * LipKRelaxation;
par.xTrue = xTrue;
par.MaxIter = MaxIter;
par.TolX = eps;
par.bSilent = false;
par.bPrimalObjectiveValue = true;
par.fhPrimalObjectiveValue = fhEnergy;
par.bCalculateRelativeError = true;
par.bPlot = false;
par.OutputInterval = 1;
par.bRelativeError = true;
par.fhRelativeError = @funRelativeL2Error;

%% BOS
% ===============================
% BOS using Xiaojing's code
% ===============================
sDescription = 'BOS';
if any(strcmp(sDescription, lAlg))
    beta = 10;
    wTV = wReg;
    wP = wTV*beta;

    maxiter = MaxIter;
    relchg_tol = eps;
    bprint = true;

    opts = [];

    opts.u0 = u0; 
    opts.bprint = bprint;
    opts.fhEnergy = fhEnergy;

    disp('TVL2: BOS algorithm is running...');
    tic;
    [ubos, outbos] = modTVL2_BOS(f, A, m, n, wTV, wP, 1, maxiter, relchg_tol, opts);
    toc;
    xBOS = ubos;
    etcBOS = [];
    etcBOS.RelativeError = outbos.err;
    etcBOS.PrimalObjectiveValue = outbos.obj;
    etcBOS.CPUTime = outbos.cpu;
end

%% BOSVS
% ===============================
% BOSVS using Xiaojing's code
% ===============================
sDescription = 'BOSVS';
if any(strcmp(sDescription, lAlg))
    beta = 10;
    wTV = wReg;
    wP = wTV*beta;

    delta = 1;

    maxiter = MaxIter;
    relchg_tol = eps;
    bprint = true;

    eta = 3;
    sigma = 0.9999;
    tau = 2;
    bb_min = 1e-3;

    opts = [];
    opts.u0 = xTrue; 
    opts.bprint = bprint;
    opts.BOS = false;
    opts.fhEnergy = fhEnergy;

    fprintf('Data undersampling ratio: %4.3g.\n', sum(p(:))/numel(p));
    disp('BOSVS algorithm is solving TVL2 problem ...');
    tic
    [ubosvs, outbosvs] = modTVL2_BOSVS(b, A, nRow, nCol, wTV, wP, delta, maxiter, relchg_tol, opts, eta, sigma, tau, bb_min);
    toc;
    xBOSVS = ubosvs;
    etcBOSVS = [];
    etcBOSVS.RelativeError = outbosvs.err;
    etcBOSVS.PrimalObjectiveValue = outbosvs.obj;
    etcBOSVS.CPUTime = outbosvs.cpu;
end

%% USL
% ===============================
% USL using Wei's code
% ===============================
sDescription = 'USL';
if any(strcmp(sDescription, lAlg))
    data = [];
    data.u0 = xTrue;
    [data.pm, data.pn]=size(data.u0);
    data.n= data.pm * data.pn;
    data.m= nSample;
    data.R = 150.0;
    data.b = f;
    data.sp = sense_map;
    data.norm = norm(xTrue, 'fro');
    data.eta = 0;
    control.theta = 0.6;
    control.beta = 0.6;
    control.lambda=wReg;
    control.epsilon = 1e-7;
    control.iter_limit = MaxIter;  %%% choose from 20 to 500.
    control.bundle_limit = 5;  %%% choose from 2 to 10.
    control.M2 = perm_gen (control.bundle_limit+2);
    
    [xout, etc] = funUSL_MFTV(data, control);
    xUSL = xout;
    etcUSL = etc;
    etcUSL.PrimalObjectiveValue = etc.ObjVal;
end

%% Accelerated Primal Dual (APD)
sDescription = 'APD';
sFunctionName = 'funAPD';
if any(strcmp(sDescription, lAlg))
    par.StepsizePolicy = 1;

    eval(['par', sDescription, ' = par;']);
    fprintf('TVL2: %s algorithm is running...\n', sDescription);
    [x,~,  etc] = feval(sFunctionName, fhGradG, fhK, fhKt, par);
    eval(['x', sDescription, '=x;']);
    eval(['etc', sDescription, '=etc;']);
end

%% Linearized Primal Dual (LPD)
sDescription = 'LPD';
sFunctionName = 'funAPD';
if any(strcmp(sDescription, lAlg))
    par.StepsizePolicy = 0;

    eval(['par', sDescription, ' = par;']);
    fprintf('TVL2: %s algorithm is running...\n', sDescription);
    [x, ~, etc] = feval(sFunctionName, fhGradG, fhK, fhKt, par);
    eval(['x', sDescription, '=x;']);
    eval(['etc', sDescription, '=etc;']);
end

%% Show results
% All available algorithms
lAlgAll = {'USL', 'LPD', 'PDx', 'NestS', 'APDx', 'APD', 'BOS', 'SBB_M', 'SBB_NM', 'ALADMM', 'ALPADMM', 'BOSVS', 'LPADMM', 'LADMM'};
lTitleAll = {'USL', 'Linearized primal-dual', 'Chambolle-Pock(x)', 'Nesterov', 'APDx', 'APD', 'BOS', 'SBB w\ linesearch', 'SBB', 'AL-ADMM', 'ALP-ADMM', 'BOSVS', 'LP-ADMM', 'L-ADMM'};

bCompareObjVal = 1;
bCompareImage = 0;
bCompareRelErr = 1;
bShowLS = 0;
% sXLabel = 'Iteration';
sXLabel = 'CPUTime';

[ism, ind] = ismember(lAlg, lAlgAll);
lAlg = lAlg(ism);
ind = ind(ism);
lTitle = lTitleAll(ind);

% xTrue = u0;
script_ComparisonPlot;