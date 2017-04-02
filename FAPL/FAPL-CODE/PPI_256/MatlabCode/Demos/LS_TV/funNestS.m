function [xagnew, etc] = funNestS(fhGradG, fhK, fhKt, par)
% =========================================================================
% =============== Write Help Text Here ====================================
% Solver for minmax 1/2 E_s[G(x, s)] + <Kx, y> - J1(y) + J2(x), or
%            minmax 1/2 G(x) + <Kx, y> - J1(y) + J2(x).

% Required input:
% fhGradG: Deterministic/stochastic gradient of G, called as fhGradG(x)
% fhK: The operator K, called as fhK(x)
% fhKt: The operator Kt, called as fhKt(x)

% Required Parameters (in par):
% LipG: Lipschitz constant for GradG
% LipK: Lipschitz constant for K
% xsize: Size of x
% ysize: Size of y

% Optional Parameters (in par):
% fhG: Function evaluation, called as fhG(x)
% fhEnergy: Energy function. If fhEnergy is not provided, the energy will
%   be calculated simply by 1/2*fhG(x) + J2(x) + J1*(Kx). The energy
%   function is used mainly for stochastic optimization purposes.
% =============== Help Text Ends ==========================================
% =========================================================================
% Notes:
% 1. The description of parameters needs to be finished.
% 2. This is a uniform template.

% --------------------------------------
% General parameters (same for all functions)
% --------------------------------------
% Required parameters
if ~exist('par', 'var')
    par = [];
end
xsize = par.xsize;
LipG = par.LipG;
LipK = par.LipK;

% % Optional parameters
% xTrue = check_par(par, 'xTrue', []);
% MaxIter = check_par(par, 'MaxIter', 100);
% bVerbose = funCheckPar(par, 'bVerbose', true);
% fhRelativeError = funCheckPar(par, 'fhRelativeError', @funRelativeL2Error);
% bRelativeError = funCheckPar(par, 'bRelativeError', false) & ~isempty(fhRelativeError) & ~isempty(xTrue);
% [bPrimalObjectiveValue, fhPrimalObjectiveValue] = funCheckPair(par, ...
%     'bPrimalObjectiveValue', 'fhPrimalObjectiveValue');
% [bDualObjectiveValue, fhDualObjectiveValue] = funCheckPair(par, ...
%     'bDualObjectiveValue', 'fhDualObjectiveValue');
% bPlot = check_par(par, 'bPlot', false);
% fhPlot = check_par(par, 'funPlot', @funPlot);
% OutputInterval = check_par(par, 'OutputInterval', 100);
% DXYRatio = check_par(par, 'DXYRatio', 1);
% StepPolicy = check_par(par, 'StepsizePolicy', 1);
% x0 = check_par(par, 'x0', zeros(xsize));
% bAverageSolution = check_par(par, 'bAverageSolution', 0);
% fhProjx = check_par(par, 'funProjx', @(x, dx)(x - dx));
% fhProjy = check_par(par, 'funProjy', @(y, dy)(y + dy));

% --------------------------------------
% Optional parameters
% --------------------------------------
% Values
MaxIter = funCheckPar(par, 'MaxIter', 100);
x0 = funCheckPar(par, 'x0', zeros(xsize));
OutputInterval = funCheckPar(par, 'OutputInterval', 1);
DXYRatio = funCheckPar(par, 'DXYRatio', 1);
StepPolicy = funCheckPar(par, 'StepsizePolicy', 1);
xTrue = funCheckPar(par, 'xTrue', []);
% Flags & Function handles
bVerbose = funCheckPar(par, 'bVerbose', true);
fhRelativeError = funCheckPar(par, 'fhRelativeError', @funRelativeL2Error);
bRelativeError = funCheckPar(par, 'bRelativeError', false) & ~isempty(fhRelativeError) & ~isempty(xTrue);
fhPlot = funCheckPar(par, 'fhPlot', @funPlot);
bPlot = funCheckPar(par, 'bPlot', false) & ~isempty(fhPlot);
fhProjx = funCheckPar(par, 'fhProjx', @(x, dx)(x - dx));
fhProjy = funCheckPar(par, 'fhProjy', @funProxMapEuclL21);
[bPrimalObjectiveValue, fhPrimalObjectiveValue] = funCheckPair(par, ...
    'bPrimalObjectiveValue', 'fhPrimalObjectiveValue');
[bDualObjectiveValue, fhDualObjectiveValue] = funCheckPair(par, ...
    'bDualObjectiveValue', 'fhDualObjectiveValue');

% --------------------------------------
% Initialization
% --------------------------------------
etc = [];
etc.CPUTime = nan(MaxIter, 1); 
etc.RelativeError = nan(MaxIter, 1); 
etc.DualityGap = nan(MaxIter, 1); 
etc.PrimalObjectiveValue = nan(MaxIter, 1); 
etc.DualObjectiveValue = nan(MaxIter, 1); 
% etc.PrimalStepsize = PrimalStepsize;
% etc.DualStepsize = DualStepsize;
% etc.AuxiliaryStepsize = AuxiliaryStepsize;
xnew = x0;
xagnew = x0;

tStart = tic;

% Stepsize
tlist = 1 : MaxIter;
switch StepPolicy
    case 1
        % Only depend on t
        alpha = 2 ./ (tlist + 1);
        mu = sqrt(1) * LipK / (MaxIter + 1) * DXYRatio;
        % TODO: check The coefficient for mu (2? sqrt(2)? or 1?)
        PrimalStepsize = (tlist + 1) ./ (2 * (LipG + 1./mu * LipK^2));
        DualStepsize = 1 / mu * ones(size(tlist));
    otherwise
        error('Unknown stepsize policy.');
end

% Save StepPolicy
etc.PrimalStepsize = PrimalStepsize;

for t = 1:MaxIter
    % --------------------------------------
    % Main iteration
    % --------------------------------------
    % ------Variable updating
    xag = xagnew;
    x = xnew;
    
    % ------Middle step
    xmd = (1 - alpha(t)) * xag + alpha(t) * x;

    % ------Dual iteration (or smoothing iteration)
    ynew = fhProjy(fhK(xmd) / mu, 0);

    % ------Primal iteration
    xnew = fhProjx(x, PrimalStepsize(t) * (fhKt(ynew) + fhGradG(xmd)));

    % ------Aggregate step
    xagnew = (1 - alpha(t)) * xag + alpha(t) * xnew;
    
    % --------------------------------------
    % Save CPU time
    % --------------------------------------
    etc.CPUTime(t) = toc(tStart);

    % --------------------------------------
    % Runtime outputs
    % --------------------------------------
    if mod(t, OutputInterval) == 0
        % Calculate the primal objective, dual objective and duality gap
        if bPrimalObjectiveValue
            etc.PrimalObjectiveValue(t) = fhPrimalObjectiveValue(xagnew);
        end
        if bPrimalObjectiveValue && bDualObjectiveValue
            etc.DualityGap(t) = etc.PrimalObjectiveValue(t) - etc.DualObjectiveValue(t);
        end
        % Calculate primal relative error to ground truth
        if bRelativeError
            etc.RelativeError(t) = fhRelativeError(xagnew, xTrue);
        end
        fprintf( 't=%d,POBJ=%e,DOBJ=%e,DualityGap=%e,RelErr=%e\n', ...
            t, etc.PrimalObjectiveValue(t), etc.DualObjectiveValue(t), etc.DualityGap(t), etc.RelativeError(t));
        % Plot
        if bPlot 
            parPlot.t = t;
            if isempty(xTrue)
                parPlot.Residual = abs(xagnew - xag);
                parPlot.sResidualTitle = [etcTermination.sChange, '=', num2str(etcTermination.Change)];
            else
                parPlot.Residual = abs(xagnew - xTrue);
                parPlot.sResidualTitle = sprintf('Relerr=%s', fhRelativeError(xagnew, xTrue));
            end
            hImage = fhPlot(xagnew, parPlot);
            parPlot.hImage = hImage;
        end
    end
end

% --------------------------------------
% Save iteration numbers
% --------------------------------------
etc.TotalIteration = t;

end