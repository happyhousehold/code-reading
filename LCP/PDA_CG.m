%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MATLAB Code for the Conditional Gradient (CG) algorithm            %
%     Author: Guanghui (George) Lan                                      %
%     Institute: University of Florida, Industrial & Systems Engineering %
%     @All rights reserved 2013                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this file implements an CG algorithm for solving QP %
function PDA_CG(control, domain, data, x0)

%%% initialization
start_CG=clock;
nstep = control.iter_limit;
% define the bundle structure

xt = x0; %GetInitPt(domain.type, domain.n, domain.R);
yt = xt;

LB = -1e40;
UB = 1e40;

sumconst = 0;
sumg = 0;
sumnu = 0;
frep = fopen('report.txt','a');
pre_xt = xt;
for k = 1: control.iter_limit
    %%%%%% Compute function value and gradient at the point y
    gammat = 2 /(k+1);
    zt = (1-gammat) * yt + gammat * xt;
        
    nuk = k;
    
    [fzt, gzt] = FirstOrderOracleQP(data,domain,zt);
    sumg = sumg + nuk * gzt;
    sumconst = sumconst + nuk *(fzt - gzt'*zt);
    sumnu = sumnu + nuk;
        
    [xt, val] = LMO(sumg/sumnu, domain);
    
    dist(k) = norm(xt-pre_xt)^2;
    avgdist(k)= sum(dist)/k;
    pre_xt = xt;
    
    LB = max(LB, val + sumconst/sumnu);
    
    if control.optstepsize == 0 %% not using optimal stepsize
        alphat = gammat;
    else
        alphat = OptimizeStepsize(data,domain,xt, yt);
    end
    yt = (1 - alphat) * yt + alphat * xt;
    
    [fyt, gyt] = FirstOrderOracleQP(data,domain,yt);

    UB = min(UB,min(fyt,fzt));
    
    if mod(k,10) == 0
 %        disp(sprintf('k=%d, UB=%.6e, LB=%.6e, Gap=%.6e\n',...
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
fprintf(frep,'**Summary: %s\n',str);
fclose(frep);
%avgdist




