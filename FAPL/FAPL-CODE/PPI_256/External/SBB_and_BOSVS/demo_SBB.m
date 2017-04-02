%% test SBB algorithm
% Xiaojing Ye (c). All rights reserved.
% Comments or questions? Email: xiaojing.ye.cn@gmail.com

% This software implements the algorithm developed in our papers:

% Computational Acceleration for MR Image Reconstruction in Partially Parallel Imaging
% X. Ye, Y. Chen, and F. Huang
% IEEE Transactions on Medical Imaging, 30(5), pp. 1055-1063, 2011.

% Bregman Operator Splitting With Variable Stepsize for Total Variation Image Reconstruction
% Y. Chen, W. W. Hager, M. Yashtini, X. Ye, and H. Zhang.
% Computational Optimization and Applications, 54(2), pp. 317-342, 2013.

% DISCLAIMER:  This code is for academic (non-commercial) use only.
% This code comes with absolutely NO warranty of any kind.

%% load data and toolboxes

% close all; clear;
addpath(genpath(pwd));
load data;

%% set common parameters for the algorithms
wTV = 1e-4;
wP = 10*wTV;
delta = 1;

maxiter = 100;
relchg_tol = 1e-5;

%% set backtracking parameters
eta = 3;
sigma = 0.9999;
tau = 2;
bb_min = 1e-3;

%% set data and operators
[m, n, k] = size(sense_map);

% set operators
AO = @(x) A_oper(x, p, m, n, k, sense_map);
AOT = @(x) At_oper(x, p, m, n, k, sense_map);
A = A_operator(@(x) AO(x), @(x) AOT(x));

% back projection image
bpr = zeros(m,n);
for ch=1:k, bpr = bpr + abs(ifft2(f(:,:,ch))).^2; end; bpr = sqrt(m*n*bpr);

%% set optional inputs
opts = [];
opts.bpr = bpr;
opts.BOS = false;

%% run SBB
fprintf('Data undersampling ratio: %4.3g.\n', sum(p(:))/numel(p));
disp('SBB algorithm is solving TVL2 problem ...');
tic
[u1, out1] = SBB(f, A, m, n, wTV, wP, delta, maxiter, relchg_tol, opts);
toc

%% run BOSVS
fprintf('Data undersampling ratio: %4.3g.\n', sum(p(:))/numel(p));
disp('BOSVS algorithm is solving TVL2 problem ...');
tic
[u2, out2] = BOSVS(f, A, m, n, wTV, wP, delta, maxiter, relchg_tol, opts, eta, sigma, tau, bb_min);
toc

%% show results

imshow(abs([bpr,u1,u2])); 
title('left: back projection, middle: SBB reconstruction, right: BOSVS reconstruction');

%% show comparison
figure;
semilogy(out1.cpu,out1.obj,out2.cpu,out2.obj);
xlabel('CPU time (s)'); ylabel('Objective');
legend('SBB (u_1)','BOSVS (u_2)');