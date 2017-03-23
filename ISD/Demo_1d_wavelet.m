% A simple demo of iterative support detection (ISD) on wavelet-sparse signal recovery

% Purpose: try to recover true signal xs from the following observations:
%   b=A*ys, where A is the measurement matrix, ys is the true signal to be 
%   recovered.
%
% ys is sparse in the Wavelet (Haar) domain. xs = W(ys) is sparse, where W is 1D Haar wavelet transform.
% Please download the required Rice Wavelet Toolbox from http://www-dsp.rice.edu/software/rice-wavelet-toolbox.

% Written by Yilun Wang and Wotao Yin 2009. Rice University

close all; clear; clc;

seed = round(5000*rand);
fprintf('Seed = %d\n',seed);
if exist('RandStream','file')
    RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed));
else
    rand('state',seed); randn('state',seed^2);
end

% Size of problems
n = 200; % Length of signals
m = 60;  % Number of measurements

% Number of nonzerom components.
% Number of nonzerom components.
switch 'difficult'
    case 'difficult'
        % DIFFICULT
        k = 22;
        disp('Test a difficult compressive sensing problem.');
    case 'easy'    
        % EASY
        k = 10;
        disp('Test an easier compressive sensing problem.');
end

% Get 1D Haar Wavelets basis using Rice Wavelets Toolbox
if ~(exist('daubcqf','file') && exist('mdwt','file') && exist('midwt','file'))
    error('Please add Rice Wavelet Toolbox to path');
end
h2 = daubcqf(2);
W.times = @(x) mdwt(x,h2);
W.trans = @(y) midwt(y,h2);

opts.basis = W;

% Generate true sparse signals
r2=randperm(n); r2=r2(1:k);
xs=zeros(n,1); xs(r2)=randn(k,1); % generate sparse Wavelets coefficients
ys=W.trans(xs);     % original signal, sparse in the Wavelets domain

%  Gaussian measurement matrix
A=randn(m,n);

%  oberservation
b=A*ys;

%  Stoping tolerance: the relative error of recovery in terms of 2-norm is
%  smaller than tol, then stop the algorithm.
opts.tol=100*eps;

%  Maximal outside iteration number:
opts.maxit=10;

%  Send the true sparse solution merely for information display
opts.xs = xs;

%  Turn debug information display on (also pauses each iteration)
% opts.Dflag = true;

%  Call ISD_Threshold
[x,Out] = Threshold_ISD_1D(A,b,opts);

% fprintf('-----------------------\n'); 
% [x,Out] = Truncated_L1_ISD(A,b,opts); % An older implementation

% Reconstruction errors in 2 norm
relerr2 = norm(x-xs,2)/norm(xs,2);
fprintf('Output relative error in 2-norm: %4.2e\n', relerr2);

if relerr2 < 1e-4, 
      fprintf('ISD: Success after solving 1 L1 and %d Truncated-L1 Min Problems\n', Out.iter-1);
else
      disp('ISD Threshold: failed. Please use less agressive threshold.');
end
