function [u, out] = SBB(f, A, m, n, wTV, wP, delta, maxiter, relchg_tol, opts)
% SBB algorithm for TVL2 image reconstruction:
% minimize_u (wTV)*TV(u) + 0.5*|Au-f|^2
%
% Input:
%   f...        data
%   A...        A operator (change definition in main test function)
%   wTV...      weight of TV regularization
%   wP...       weight of penalty term
%   opts...     other parameter options
%
% out:
%   u...        reconstruction
%   out...   outs, e.g. objective value

%% check input
if (wTV <= 0); error('Weight parameter of TV term must be positive'); end;

%% initialize variables

% pre-compute F'(D'D)F
rho = wTV*wP; sqrt_rho = sqrt(rho);
denom = abs(psf2otf([sqrt_rho,-sqrt_rho],[m,n])).^2 + ...
    abs(psf2otf([sqrt_rho;-sqrt_rho],[m,n])).^2;

% claim variable spaces, use 0 as initial guess
u = zeros(m,n); [ux, uy] = grad2(u);
wx = zeros(m,n); wy = zeros(m,n);
bx = zeros(m,n); by = zeros(m,n);
r = -f; % residual r = Au - f
At_r = A'*r; % conjugate of A applied to r

% outs
out = []; out.obj = []; out.err = []; out.cpu = [];

%% main loop
t0 = cputime;
for iter = 1:maxiter
    
    % -------------------------------------------------
    % iteration log
    %
    out.cpu = [out.cpu; cputime-t0];
    out.obj = [out.obj; objval(u, ux, uy, r, wTV)];
    if isfield(opts,'u0'), out.err = [out.err; relerr(u, opts.u0)]; end;
        
    % -------------------------------------------------
    % u - subproblem
    %
    Dt_wb = compute_Dt_wb(wx,wy,bx,by,rho);
    u_prev = u; 
    At_r_prev = At_r;
    u = ifft2(fft2(Dt_wb + delta*u_prev - At_r_prev)./ (denom+delta));
    
    [ux, uy] = grad2(u);
    r = A*u - f;
    At_r = A'*r;

    % -------------------------------------------------
    % w-subproblem
    %
    [wx, wy] = compute_w(ux,uy,bx,by,rho,wP);
        
    % -------------------------------------------------
    % check stopping criterion
    %
    relchg = norm(u - u_prev,'fro')/norm(u,'fro');
    if relchg < relchg_tol
        break;
    end
    
    % -------------------------------------------------
    % update delta
    %
    chg_u = u - u_prev;
    chg_At_r = At_r - At_r_prev;
        
    bbnumer = sum(sum(conj(chg_At_r).*chg_u));
    delta = real(bbnumer/norm(chg_u,'fro')^2);

    if (delta < 1e-10) || (delta > 1e+10)
        fprintf('Warning: BB stepsize is too small/big: %e\n', delta)
        delta = 1;
    end
    
    if opts.BOS, delta=0.5; end;
    
    % -------------------------------------------------
    % multiplier update
    %
    bx = bx + rho*(ux - wx);
    by = by + rho*(uy - wy);
    
end

out.iter = iter;