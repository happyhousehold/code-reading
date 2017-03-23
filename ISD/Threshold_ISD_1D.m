function  [x,Out] = Threshold_ISD_1D(A,b,opts)
% Thesholding-based iterative support detection (ISD) for 1D signals


%%%%%%%%%% Input specification %%%%%%%%%%%%%%%%%%%
% A: measurement matrix or operator denoted by a structure.
% b: measured data (could be noisy)
% opts (optional):
%     opts.xs:      true signal; if provided, debug flag (Dflag) can be set to true
%     opts.sigma:   standard variation of added zero mean Gaussian noise to b.
%
%     opts.rho:     parameter for YALL1: regularization parameter for min |x|_1 + 1/2/rho \|Ax-b\|^2
%                   see YALL1's manual.
%     opts.tol:     final stopping parameter for YALL1; see YALL1's manual.
%    opts.loose_tol:intermediate stopping parameter for YALL1; see YALL1's manual.
%
%     opts.basis:   sparsifying basis such as wavelets. By default, it is the identity matrix if not provided.
%                   Assume the basis is W. Then opts.basis must have two fields:
%                            opts.basis.times: a function handle for W*x
%                            opts.basis.trans: a function handle for W'*x
%                   This option is used by both ISD's main code and YALL1. See YALL1's manual.
%
%     opts.suppDetct:   1 - ISD (default); 0 - (NOT YET TESTED) reweighted L1 minimization by Cades et al.
%     opts.alpha:   parameter for support detection; see nested function compute_threshold() below
%     opts.maxit:   maximal number of thresholding (or reweighted) iterations. 9 by default
%                   if 1, then just plain L1 minimization.
%
%     opts.Dflag:           flag about whether diagnostic plots of output of each innter iteration will be drawn
%
%     opts.signalSize:      size of true signal;
%
%%%%%%%%%% Input specification %%%%%%%%%%%%%%%%
%
% x: reconstructed signal
% Out: (optional)
%     Out.iter: total number of ISD iterations
% %     Out.d: if Out.d(itr)==1, successful reconstructoin on the itr iteration. 
% if opts.xs is provided, then the following fields may also exist
%     Out.SNR:           iterative SNRs
%     Out.Error_inf:     absolute error in inf-norm. 
%     Out.RelError2:     relative error in l2 norm. 
%     Out.RelError1:     relative error in l1 norm. 

% Written by Yilun Wang and Wotao Yin Nov 2009
%
% YALL1: Copyright by Yin Zhang, http://www.caam.rice.edu/~optimization/L1/YALL1/


%% Option loading and variable initialization
if ~exist('opts','var'); opts = []; end

% ISD options
if isfield(opts,'maxit'); maxit = opts.maxit; else maxit = 10; end
if isfield(opts,'xs');
    xs=opts.xs(:);
    SNR=zeros(maxit,1); snr([],xs);
    err_inf=zeros(maxit,1);
    relerr2=zeros(maxit,1);
    relerr1=zeros(maxit,1);
    nrm_xs = norm(xs); nrm1_xs = norm(xs,1);
end
if isfield(opts,'sigma'); sigma = opts.sigma; else sigma = 0; end
if isfield(opts,'alpha'); alpha = opts.alpha; else alpha = 7; end
if isfield(opts,'suppDetct');   suppDetct = opts.suppDetct; else suppDetct = 1; end
if isfield(opts,'Dflag') && opts.Dflag==1 && exist('xs','var') && suppDetct;  Dflag =1; else Dflag = 0; end

% YALL1 options
if isfield(opts,'rho');   rho = opts.rho; else  rho = 0; end
if isfield(opts,'basis'); basis = opts.basis; end
if isfield(opts,'tol');   tol = opts.tol; else tol = 1e-6; end
if isfield(opts,'loose_tol'); loose_tol = opts.loose_tol; else loose_tol = 1e-2; end
% YALL1 input/output initialization
opts1=[]; yall1_out=[];

% general initialization
if ~isstruct(A); n=size(A,2); else  n=A.n; end

m=length(b); nrm_b = norm(b);
W=ones(n,1);    % initial weights
x=zeros(n,1);
thresh=inf;     % support detection threshold

bBreak = 0;   % 0 - continue; 1 - stop after another L1 min; 2 - stop immediately
isemp  = 0;   % Whether two sequenct detected supports are indentical; 1 - last two are equal, 2 - last three are equal
 
 
%% Main iteration part
for itr=1:maxit
    
    % Solve truncated L1 minimization
    Solve_Truncated_L1(); % generate a new sparse point x
    
    if bBreak == 1 || itr == maxit; 
        if exist('xs','var'); Inner_record(false); end
        break; 
    end

    % Support detection
    W_prev = W;
    bBreak = support_detection(); % updates weights and check stopping rules
    if bBreak == 2; break; end
    
    % If xs is provided, compute and record the errors
    if exist('xs','var'); Inner_record(); end
end

%% Post processing
Out.iter=itr;
if exist('xs','var'); Final_record_RelErr(); end
 
%% Nested functions

% Record relative errors
    function Inner_record(bShowDet)
        
        SNR(itr)=snr(x);
        err_inf(itr) = norm(x-xs,'inf');
        relerr2(itr) = norm(x-xs)/nrm_xs;
        relerr1(itr) = norm(x-xs,1)/nrm1_xs;
        
        fprintf('itr=%d inf=%4.2e 2norm=%4.2e 1norm=%4.2e ', itr, err_inf(itr), relerr2(itr), relerr1(itr));
        
        nozero_criteria = 1e-6;
        
        true_nzs = (abs(xs)>nozero_criteria);
        nzs = (abs(x)>nozero_criteria); % jumps is logical
        good_nzs = nzs & true_nzs;
        bad_nzs = nzs & (~true_nzs);
        miss_nzs = true_nzs & ~nzs;

        fprintf('x_supp(total,good,bad,miss)=(%2d,%2d,%2d,%2d) ',nnz(nzs),nnz(good_nzs),nnz(bad_nzs),nnz(miss_nzs));

        if nargin==0 || bShowDet
            det_nzs = ~W;
            det_good_nzs = det_nzs & true_nzs;
            det_bad_nzs = det_nzs & (~true_nzs);
            
            fprintf('det_supp(total,good,bad)=(%2d,%2d,%2d) ',nnz(det_nzs),nnz(det_good_nzs),nnz(det_bad_nzs));
        end
        
        fprintf('\n');
        
        % Show debug information
        if Dflag;
            figure(itr);

                stem(1:n,xs,'k.','MarkerSize',15); hold on;
                stem(find(good_nzs),x(good_nzs),'ro','LineWidth',2);
                if any(bad_nzs)
                    stem(find(bad_nzs), x(bad_nzs),'bd','LineWidth',2);
                    legend('true signal','true nonzero','false nonzero');
                else
                    legend('true signal','true nonzero');
                end
                plot([1 n],[thresh thresh],'--g');
                plot([1 n],[-thresh -thresh],'--g');
                hold off
                title(sprintf('Itr=%d (total,good,bad,miss)=(%2d,%2d,%2d,%2d) RelErr=%4.2e',...
                      itr, nnz(nzs),nnz(good_nzs),nnz(bad_nzs),nnz(miss_nzs),relerr2(itr)),...
                      'fontsize',14);
              
            disp('Paused. Press any key to continue...');
            pause;
        end

    end

    function Final_record_RelErr()
        Out.SNR = SNR;
        Out.Error_inf = err_inf;
        Out.RelError2 = relerr2;
        Out.RelError1 = relerr1;
    end

% Solve truncated L1 minimization
    function Solve_Truncated_L1()
        
        opts1.rho=rho;
        opts1.tol=tol;
        opts1.maxit=1500;
        if suppDetct && maxit>1
            switch itr
                case 1;     
                    opts1.tol=1e-1;
                case maxit;
                    opts1.tol=tol;
                    opts1.maxit=3000;
                otherwise;
                    if bBreak == 1;
                        opts1.maxit=3000;
                        opts1.tol=tol;
                    else
                        opts1.tol=loose_tol;
                    end
            end
        end
        opts1.print = 0;
        opts1.weights = W;
        if exist('basis','var');   opts1.basis = basis;  end

        [x, yall1_out] = yall1_ext(A,b,opts1);
        
        % for yall1 warm start
        opts1.x0 = x;
        opts1.bmax0 = yall1_out.bmax;
        opts1.y0 = yall1_out.y;
        opts1.z0 = yall1_out.z;
        
    end

% Nested fuction for support detection for ISD or computing weights for iterative
% reweighted L1 algorithm
    function res = support_detection() % For both IRL1 and ISD
        %  bBreak: 0 - continue; 1 - stop after another L1 min; 2 - stop immediately
        
        res = 0;

        switch suppDetct
            case 0 % reweighting presented in  "Enhancing Sparsity by Reweighted l1 Minimization"
                % by E. J. Candes, M. B. Wakin, and S. Boyd.
                min_RegPara=max(sigma, 2^(-8));
                W=1./(abs(x)+max(1/2^(itr-1), min_RegPara));
                W=W/max(W);
                
            case 1 % Truncated L1 minimization---First Significant Jump Method.
                
                x_abs = abs(x);
                x2_abs = sort(x_abs); % x2_abs = x2_abs(x2_abs>eps);
                
                % if the tail is small enough, STOP
                if (sigma == 0 && x2_abs(n-round(m/2))/nrm_b < 1e-14) || ...  % noiseless case
                   (sigma > 0  && x2_abs(n-round(m/2))/nrm_b < max(tol,sigma/100));% noisy case
                    W = []; res = 1;
                end
                
                % compute threshold based on the first significant jump rule
                difference=diff(x2_abs);
                xmax=x2_abs(end);
                thresh = compute_threshold(alpha);
                
                % compute weights
                W=ones(size(x));
                W(x_abs>thresh)=0;
                
                % if W contains too few nonzeros, STOP after another L1 minimization
                if nnz(W) < n-m; 
                    res = 1; return; 
                end
                
                W_diff_tol = (n-norm(W_prev,1))/10;
                if norm(W-W_prev, 1) < W_diff_tol
                    % We expect at least 10% update in W.
                    % However, W may be almost equals W_prev due to
                    % (1) solution is already sparse. reducing threshold won't make a difference
                    % (2) the loose stopping tolerance of YALL1.
                    % (3) threshold is too loose.
                    % SOLUTION: 
                    % first time, reduce YALL1's tolerance
                    % second time, reduce alpha to see if a different W can be obtained; if not, force stop after another L1
                    % third time, force stop after another L1
                    isemp = isemp + 1;
                    switch isemp
                        case 1
                            loose_tol = 0.01*loose_tol;
                        case 2
                            thresh = compute_threshold(alpha,ceil(W_diff_tol));
                            W=ones(size(x));
                            W(x_abs>=thresh)=0;
                            res = 0;
                        otherwise
                            res = 1;
                    end
                    % recompute threshold with a smaller alpha
                else
                    isemp = 0; res = 0;
                end
                
            otherwise
                error('ISD_method has wrong value');
        end
        
        function thresh1 = compute_threshold(a,bMore)
            
            difftol = a*xmax/m/min(itr,6); % sparse Gaussian signals.
            %difftol =  8*xmax/m/20/itr; % sparse signles. pr=1/3
            %difftol =  8*xmax/m/40/itr; % compressible signles. pr=1/3
            %difftol =  8*xmax/m/5/itr;  % sparse signles. pr=0.8
            %difftol =  8*xmax/m/5/itr;  % sparse signles. pr=1
            %difftol =  8*xmax/m/2/itr;  % sparse signals. pr=2
            %difftol =  8*xmax/m/2/itr;  % sparse signals. pr=4
            %difftol =  8*xmax/m/itr;    % Phantom images:
            %difftol =  8*xmax/m/itr;      or other options ? % Wavelets of Cameraman:
            %difftol =  8*xmax/m/itr;    % Phantom images: xmeth=5;
            
            large = difference>difftol;
            
            if any(large)
                pos = find(large,1,'first');
                if nargin == 2; pos = max(pos - bMore,1); end
                thresh1=x2_abs(pos+1);
            else
                thresh1=x2_abs(1);
            end
        end
        
    end


end % main function






    
     
