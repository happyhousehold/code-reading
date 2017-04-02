function [xag, etc] = funLS_TV_APDU(A, b, wTV, par)
% =========================================================================
% =============== Write Help Text Here ====================================
% =============== Help Text Ends ==========================================
% =========================================================================
% Notes:
% 1. The description of parameters needs to be finished.

LipG = par.LipG;
LipK = par.LipK;
xsize = par.xsize;
xTrue = par.xTrue;
bVerbose = funCheckPar(par, 'bVerbose', true);
chi = double(funCheckPar(par, 'bPreconditioned', 0));
MaxIter = funCheckPar(par, 'MaxIter', 100);
RhoEst = funCheckPar(par, 'RhoEstimate', 1/LipK);

FKtKFt = wTV^2 .* (abs(psf2otf([1,-1],xsize)).^2 + abs(psf2otf([1;-1],xsize)).^2);

% --------------------------------------
% Initialization
% --------------------------------------
etc = [];
etc.CPUTime = nan(MaxIter, 1);
etc.RelativeError = nan(MaxIter, 1);
etc.PrimalObjectiveValue = nan(MaxIter, 1);

x = zeros(xsize);
xag = zeros(xsize);
w = zeros([xsize, 2]);
y = zeros([xsize, 2]);
Kx = zeros([xsize, 2]);
tStart = tic;

for t = 1:MaxIter
    % --------------------------------------
    % Main iteration
    % --------------------------------------
    % ------Stepsize parameters
    alpha = 2/(t+1);
    rho = t/MaxIter*RhoEst;
    tau = MaxIter/t*RhoEst;
    theta = tau;
    eta = (2*LipG + chi*MaxIter*LipK^2)/t;
    
    % ------Middle step
    xmd = (1 - alpha) * xag + alpha * x;
    gradx = A'*(A*xmd - b);
    
    % ------x iteration (with or without preconditioning)
    if chi
        x = x - (funTVNegDiv(theta*(Kx - w) + y, wTV, 1) + gradx)/eta;
    else
        %   x_{t+1} = argmin theta_t*|Kx - Bw_t + b|^2/2 + <gradG, x> + <y_t, Kx> + eta_t*|x -x_t|^2/2
        %           = argmin theta_t*|Kx|^2/2 + <x, K'(theta_t*(b-Bw_t)+y_t) + gradG> + + eta_t*|x -x_t|^2/2
        RHS = funTVNegDiv(theta*w - y, wTV, 1) - gradx + eta*x;
        LHS = theta*FKtKFt + eta;
        x = ifft2(fft2(RHS) ./ LHS);
    end
    Kx = funTVGrad(x, wTV, 1);
    % ------w iteration
    w = funThres21Norm(Kx + y/tau, tau);
    % ------y iteration
    y = y - rho * (w - Kx);
    % ------Aggregate step
    xag = (1 - alpha) * xag + alpha * x;

    % --------------------------------------
    % Save CPU time
    % --------------------------------------
    etc.CPUTime(t) = toc(tStart);
    % --------------------------------------
    % Runtime outputs
    % --------------------------------------
    % Calculate the primal objective
    etc.PrimalObjectiveValue(t) = par.fhPrimalObjectiveValue(xag);
    % Calculate primal relative error to ground truth
    etc.RelativeError(t) = funRelativeL2Error(xag, xTrue);
    funPrintf(bVerbose, 't=%d,POBJ=%e,RelErr=%e\n', ...
        t, etc.PrimalObjectiveValue(t), etc.RelativeError(t));
end

end