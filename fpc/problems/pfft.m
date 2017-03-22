function y = pfft(x,mode,n,picks)

% Implements partial fourier transform matrix defined by taking only the 
% rows in the m-vector picks.  Resulting matrix A is m x n, m < n.  pfft
% returns the following quanities based on the value of mode:
%
%  -1: y = picks
%   1: y = A*x
%   2: y = A'*x

switch mode
    case -1
        y = picks;
    case 1
        % y = A*x
        y = fft(x);
        y = y(picks)/sqrt(n);
    case 2
        % y = A'*x
        y = zeros(n,1);
        y(picks) = x;
        y = ifft(y)*sqrt(n);
    otherwise
        error('mode must be 1 (for A*x) or 2 (for A''*x).');
end

return

% Last modified 11 January 2008.