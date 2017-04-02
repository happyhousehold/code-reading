% Solver for the penalized inaccurate matrix game
%   min_x max_y .5*|Ax-b|^2 + <Kx, y>
% where x and y are on simplices, and Kx and K'y are calculated through
% calls of deterministic/stochastic oracles

% References: 
% [1] Nemirovski, A. (2004). Prox-method with rate of convergence O
% (1/t) for variational inequalities with Lipschitz continuous monotone
% operators and smooth convex-concave saddle point problems. SIAM Journal
% on Optimization, 15(1), 229-251.
% 
% [2] Juditsky, A., Nemirovskii, A. S., & Tauvel, C. (2008). Solving
% variational inequalities with stochastic mirror-prox algorithm. arXiv
% preprint arXiv:0809.0815.
function [xav, yav, etc] = funG_PM(Q, fhK, fhKt, K, par)

[yLength, xLength] = size(K);
if nargin<5
    par = [];
end
bSilent = check_par(par, 'bSilent', false);
fhDualityGap = check_par(par, 'fhDualityGap', []);
bDualityGap = check_par(par, 'bDualityGap', false) & ~isempty(fhDualityGap);
OutputInterval = check_par(par, 'OutputInterval', 1);
MaxIter = check_par(par, 'MaxIter', 100);
TolGap = check_par(par, 'TolGap', 1e-3);
M = check_par(par, 'M', 0);
LipG = check_par(par, 'LipG', 0);
LipK = check_par(par, 'LipK', 1);
theta = check_par(par, 'theta', 5);

% --------------------------------------
% Initialization
% --------------------------------------
etc = [];
etc.CPUTime = nan(MaxIter, 1); 
etc.PrimalObjectiveValue = nan(MaxIter, 1); 
etc.DualObjectiveValue = nan(MaxIter, 1); 
etc.DualityGap = nan(MaxIter, 1); 
xnew = ones(xLength, 1) / xLength;
ynew = ones(yLength, 1) / yLength;
xav = zeros(xLength, 1);
yav = zeros(yLength, 1);
L = sqrt(2*log(xLength)*(LipG + LipK)^2 + 2*log(yLength)*LipK^2);
if M ==0
    Stepsize = 1/L/sqrt(2);
else
    Stepsize = min(1/(sqrt(3)*L), 2/M/sqrt(21*MaxIter)) * theta;
end
xStep = Stepsize * 2 * log(xLength);
yStep = Stepsize * 2 * log(yLength);

tStart = tic;
for t = 1:MaxIter
    % --------------------------------------
    % Main iteration
    % --------------------------------------
    % ------Variable updating
    y = ynew;
    x = xnew;

    % ------Extragradient step    
    xeg = funProxMapEntropy(x, xStep * (Q*x + fhKt(y)));
    yeg = funProxMapEntropy(y, -yStep* fhK(x));
    
    % ------Gradient step
    xnew = funProxMapEntropy(x, xStep * (Q*x + fhKt(yeg)));
    ynew = funProxMapEntropy(y, -yStep * fhK(xeg));

    % ------Aggregate step
    xav = (xav*(t-1) + xeg)/t;
    yav = (yav*(t-1) + yeg)/t;

    % --------------------------------------
    % Save CPU time
    % --------------------------------------
    etc.CPUTime(t) = toc(tStart);
    % --------------------------------------
    % Calculate the duality gap
    % --------------------------------------
    if bDualityGap && mod(t, OutputInterval) == 0
        [etc.DualityGap(t), etc.PrimalObjectiveValue(t), etc.DualObjectiveValue(t)]...
            = fhDualityGap(xav, yav);
        silent_fprintf(bSilent, 't=%d, POBJ=%e, DOBJ=%e, DualityGap=%e\n', ...
            t, etc.PrimalObjectiveValue(t), etc.DualObjectiveValue(t), etc.DualityGap(t));
        if etc.DualityGap(t) < TolGap
            break;
        end
    end
    
end
    
% --------------------------------------
% Save total iteration number
% --------------------------------------
etc.TotalIteration = t;

end