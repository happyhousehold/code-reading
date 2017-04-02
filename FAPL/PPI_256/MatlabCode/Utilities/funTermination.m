function [bTermination, etc] = funTermination(varargin)
% Output
%   bTermination: true if it is good to stop.
%   etc.Change: Difference between xnew & x
%   etc.sChange: Describe etc.Change (Relative change, absolute change,
%      ...)

    if nargin == 2
        SVars = varargin{1};
        par = varargin{2};
        if ~exist('par', 'var')
            par = [];
        end
        sCriterion = check_par(par, 'sCriterion', 'L2RelativeChange');
        bTermination = false;
        switch sCriterion
            case 'L2RelativeChange'
                numer = norm(SVars.xsolnew(:) - SVars.xsol(:));
                denom = norm(SVars.xsol(:));
                relchg = numer / denom;
                if denom == 0 && numer == 0
                    bTermination = true;
                end
                if denom ~= 0 && relchg < SVars.TolX
                    bTermination = true;
                end
                etc.sChange = sprintf('Rel. change=%1.2e', relchg);
            otherwise
                error('Unknown stopping criterion: %s', sCriterion);
        end
    else
        xnew = varargin{1};
        x = varargin{2};
        TolX = varargin{3};
        if nargin == 4
            sCriterion = varargin{4};
        end
        if ~exist('sCriterion', 'var')
            sCriterion = 'L2RelativeChange';
        end
        bTermination = false;
        switch sCriterion
            case 'L2RelativeChange'
                numer = norm(xnew(:) - x(:));
                denom = norm(x(:));
                relchg = numer / denom;
                if denom == 0 && numer == 0
                    bTermination = true;
                end
                if denom ~= 0 && relchg < TolX
                    bTermination = true;
                end
                etc.Change = relchg;
                etc.sChange = 'Relative change';
            otherwise
                error('Unknown stopping criterion: %s', sCriterion);
        end
    end
        

end
