%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MATLAB Code for the Conditional Gradient (CG) algorithm            %
%     Author: Guanghui (George) Lan                                      %
%     Institute: University of Florida, Industrial & Systems Engineering %
%     @All rights reserved 2013                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this file implements an CG algorithm for solving QP %
function CG(control, domain, data, x0)

%%% initialization
start_CG=clock;
nstep = control.iter_limit;
% define the bundle structure

xt = x0; %GetInitPt(domain.type, domain.n, domain.R);
yt = xt;

LB = -1e40;
UB = 1e40;
frep = fopen('report.txt','a');
for k = 1: control.iter_limit
    %%%%%% Compute function value and gradient at the point y
    [fyt, gyt] = FirstOrderOracleQP(data,domain,yt);
    [xt, val] = LMO(gyt, domain);
    
    dist(k) = norm(xt-yt)^2;
    avgdist(k)= sum(dist)/k;
    
    LB = max(LB, fyt + val - gyt' * yt);
    %UB = min(UB,fyt);
    if UB > fyt
        UB = fyt;
        outx = yt;
    end
    
    if control.optstepsize == 0 %% not using optimal stepsize
        alphat = 2/(k+1);
    else
        alphat = OptimizeStepsize(data,domain,xt, yt);
    end
    yt = (1 - alphat) * yt + alphat * xt;
    
    if mod(k,10) == 0
 %       disp(sprintf('k=%d, UB=%.6e, LB=%.6e, Gap=%.6e\n',...
 %         k, full(UB), full(LB), full(UB-LB)));
        str = sprintf('k=%d, UB=%.6e, LB=%.6e, Gap=%.6e\n',...
          k, full(UB), full(LB), full(UB-LB));
        disp(str);
        fprintf(frep,str);
    end
end

time_CG=etime(clock,start_CG);
disp(sprintf('end of execution.\n'));
str  = sprintf('nstep=%d, UB=%.2e, LB=%.2e,time=%5.2f\n',...
    nstep, full(UB), full(LB), time_CG);
disp(str);
%if domain.type == 4
%    disp(sprintf('rank=', rank(reshape(outx,domain.n,domain.n),1e-6)));
%end
%frep = fopen('report.txt','a');
fprintf(frep,'**Summary: %s\n',str);
fclose(frep);
%avgdist




