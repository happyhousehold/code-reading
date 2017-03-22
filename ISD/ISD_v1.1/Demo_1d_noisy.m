% A simple demo of iterative support detection (ISD) on noisy measurements

% Purpose: try to recover true signal xs from the following observations:
%   b=A*xs, where A is the measurement matrix, xs is the true signal to be 
%   recovered.
%
% xs is sparse.

% Written by Yilun Wang and Wotao Yin  2009. Rice University

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
k = 10;

% Generate true sparse signals
r2=randperm(n);r2=r2(1:k); xs=zeros(n,1);
xs(r2)=randn(k,1);

%  Gaussian measurement matrix
A=randn(m,n);

%  oberservation
sigma = 0.01;
b=A*xs + sigma*randn(m,1);

opts.sigma = sigma;

%  Stoping tolerance: the relative error of recovery in terms of 2-norm is
%  smaller than tol, then stop the algorithm.
opts.tol=1e-2*sigma;

%  Maximal outside iteration number:
opts.maxit=10;

%  Send the true solution merely for information display
opts.xs = xs;

%  Turn debug information display on (also pauses each iteration)
% opts.Dflag = true;

%  Call ISD_Threshold
[x,Out] = Threshold_ISD_1D(A,b,opts);
relerr = norm(x-xs,1)/norm(xs,1);
fprintf('ISD: Output relative error in 1-norm: %4.2e\n\n', relerr);

%  Call Basis Pursuit
opts.maxit=1;
[x_bp,Out_bp] = Threshold_ISD_1D(A,b,opts);
relerr_bp = norm(x_bp-xs,1)/norm(xs,1);
fprintf('BP : Output relative error in 1-norm: %4.2e\n', relerr_bp);

