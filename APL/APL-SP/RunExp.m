%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MATLAB Code for the Accelerated Prox-Level (APL) algorithm            %
%     Author: Guanghui (George) Lan                                      %
%     Institute: University of Florida, Industrial & Systems Engineering %
%     @All rights reserved 2010                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this file implements an interface for testing the APL algorithm %
clear all;
name='noname';
sfile='no';
fname='noname.ins';

flag_dat = 0;

clear data;
% control parameters
control.theta = 0.5;
control.lambda = 0.9;
control.iter_limit = 500;
control.epsilon = 1e-6;
control.bundle_limit = 30;
control.lb_mode = 0; %% 0 if minimize a linear function, 1 to minimize the
                     %% piecewise linear function.
control.DomainType = 1; %% 0 for a simplex, >=1 for a box
dataF.data.diff_eig = 0;

                     
while(1==1)
    strD=sprintf('Data: %s,\n senario: %s',name,sfile);
    imenu=menu(strD,'Load data','Run APL', 'View protocol','Quit');
    if imenu==3,
        edit report.txt
    end;
    if imenu==4,
        break;
    end;
    
    if imenu == 1,
        ls *.m
        clear A1;
        clear A1_bnd;
        clear A2;
        clear A2_bnd;
        clear V1_bnd;
        clear V2_bnd;
        clear RHS_ind;
        clear senario;

        name = input('load the data > ','s');
        %str = sprintf('%s',name);
        eval(name);
        
        sfile = input('load the scenario > ', 's');
        eval(sfile);
        
        domain.n = n1;
        data.m1 = m1;
        data.n1 = n1;
        data.m2 = m2;
        data.n2 = n2;
        data.n = domain.n;
        data.dim_xi = dim_xi;
        
%        data.Ar = exAr(:, n1+1:n1+n2);
%        data.p_rhs = xi_c_ub;
        data.mu = 0;
        x = zeros(domain.n,1);
         
        data.ns = ns;
        data.senario = senario;
        %%% define the master problem
        data.master.c = A1(1, 1:n1)';
        %%% define the constraints X for the first stage
        data.master.a = A1(2:m1+1, 1:n1); % the first row is objective coefficients
        data.master.blc = zeros(m1,1);
        data.master.buc = zeros(m1,1);
        for i = 1:m1,
            if A1_bnd(i+1, 1) == 110,
                data.master.blc(i) = -inf; % free
                data.master.buc(i) = inf;
            end
            if A1_bnd(i+1, 1) == 111, % lower bound
                data.master.blc(i) = A1_bnd(i+1, 2);
                data.master.buc(i) = inf;
            end
            if A1_bnd(i+1, 1) == 112, % upper bound
                data.master.blc(i) = -inf;
                data.master.buc(i) = A1_bnd(i+1, 3);
            end
            
            if A1_bnd(i+1, 1) == 113 || A1_bnd(i+1, 1) == 114, % double-bounded or fixed 
                data.master.blc(i) = A1_bnd(i+1, 2);
                data.master.buc(i) = A1_bnd(i+1, 3);
            end
        end
        
        %%% define the  bounds for the variables
        data.master.blx = zeros(n1,1);
        data.master.bux = zeros(n1,1);
        for i = 1:n1,       
            switch V1_bnd(i,1)
                case 110
                    data.master.blx(i) = -inf; % free
                    data.master.bux(i) = inf;
                case 111 % lower bound
                    data.master.blx(i) = V1_bnd(i, 2);
                    data.master.bux(i) = inf;
                case 112 % upper bound
                    data.master.blx(i) = -inf;
                    data.master.bux(i) = V1_bnd(i, 3);
                case {113,114} % double bound or fix
                    data.master.blx(i) = V1_bnd(i, 2);
                    data.master.bux(i) = V1_bnd(i, 3);
                otherwise
                    disp('Unknown bounds.')
            end
            
        end
    
        %%%% define the second stage problem
        %%%% we assume that the second stage problem is a min problem
        %%%% we assume that both constraints are bounded from below,
        %%%% above, or fixed (equality), the variables are either free or
        %%%% bounded below from zero
        
        %% define the dual problems
        %% define the dual objective
        data.A2 = A1(m1+1+1:m1+1+m2,:);
        data.org_dual.c = zeros(m2,1);
        data.org_dual.blx = zeros(m2,1);
        data.org_dual.bux = zeros(m2,1);
        x = zeros(n1,1);
        
        data.second_rhs = zeros(1,m2);
        
        %% the bounds for the dual variables
        for i = 1:m2,
            switch A1_bnd(m1+i+1,1)
                case 111 % lower bound
                    data.second_rhs(i) = A1_bnd(m1+i+1, 2);
                    data.org_dual.blx(i) = 0;
                    data.org_dual.bux(i) = inf;
                case 112 % upper bound
                    data.second_rhs(i) = A1_bnd(m1+i+1, 3);
                    data.org_dual.blx(i) = -inf;
                    data.org_dual.bux(i) = 0;
                case 114 % equality
                    data.second_rhs(i) = A1_bnd(m1+i+1, 2);
                    data.org_dual.blx(i) = -inf;
                    data.org_dual.bux(i) = inf;
                otherwise
                    disp('Unknown bounds.')
            end
        end
        
       %% incorporate randomness into the RHS of the second stage
        for i=1:ns,
            data.rand_second_rhs(i,:) = data.second_rhs;
            for j = 1:dim_xi  
                data.rand_second_rhs(i,RHS_ind(j)-1) = senario(i,j);
            end
        end
        
        % the objective
        data.org_dual.c = data.rand_second_rhs(1,:)' - data.A2(:, 1:n1) * x;
        
        % the constraints
        data.org_dual.a = data.A2(:, n1+1:n1+n2)';
        
        % the bounds on constraints
        for i = 1:n2,
            switch V1_bnd(n1+i,1)
                case 111 % lower bound 0
                    if V1_bnd(n1+i,2) == 0,
                        data.org_dual.blc(i) = -inf;
                        data.org_dual.buc(i) = A1(1, n1 + i);
                    else
                        disp('Unknown bounds.')
                    end
                case 110 % free
                    data.org_dual.blc(i) = A1(1, n1 + i);
                    data.org_dual.buc(i) = A1(1, n1 + i);
                otherwise
                    disp('Unknown bounds.')
            end
        end

%        term20sol;
%        data.iniSol = termsol;
%        termX = A2(2:m2+1,1:n1) * termsol;
%         data.sec_prob.c = A2(1, n1+1:n1+n2)';
%         data.sec_prob.a = A2(2:m2+1, n1+1:n1+n2);
%         data.sec_prob.blc = zeros(m2,1);
%         data.sec_prob.buc = zeros(m2,1);
%         for i = 1:m2
%             if A2_bnd(i+1,1) == 112,
%                 data.sec_prob.blc(i) = -inf;
%                 data.sec_prob.buc(i) = A2_bnd(i+1, 3)-termX(i);
%             else
%                 if A2_bnd(i+1,1) == 114,
%                     data.sec_prob.blc(i) = A2_bnd(i+1, 2)-termX(i);
%                     data.sec_prob.buc(i) = A2_bnd(i+1, 3)-termX(i);
%                 end
%             end          
%         end
%         data.sec_prob.blx = zeros(n2,1);
%         data.sec_prob.bux = inf * ones(n2,1);
%           
%         [r,od_res] = mosekopt('minimize echo(0)', data.sec_prob);
        
        %% define the smoothed dual problem
        data.sm_dual.q = data.mu * eye(data.m2) /2;
        data.sm_dual.c = -data.org_dual.c;
        data.sm_dual.a = data.org_dual.a;
        data.sm_dual.blc = data.org_dual.blc;
        data.sm_dual.buc = data.org_dual.buc;
        data.sm_dual.blx = data.org_dual.blx;
        data.sm_dual.bux = data.org_dual.bux;
         
        
        flag_dat = 1;
    end;
    
    %%%%%%%%%
    if imenu==2,
        if flag_dat==0,
            error('No data specified');
            break;
        end;
        % default control options                         
        control.theta = 0.5;
        control.lambda = 0.5;
        control.iter_limit = 500;
        control.epsilon = 1.0e-6;
        control.bundle_limit = 30;
        control.lb_mode = 0; %% 0 if minimize a linear function, 1 to minimize the
                                 %% piecewise linear function.
        
        control.SubSolver = 0;  %% 0 if using IPM for the primal subproblem
                                %% 1 if using Accelerated Bundle-level
                                %% method for the dual subproblem

        %%%%% Run 
        while(1==1)
            
            % choose options
            strP = sprintf('Processing data %s',name);
            strIter = sprintf('iteration limit: %d', control.iter_limit);
            strEpsilon = sprintf('epsilon: %.2e', control.epsilon);
            strTheta = sprintf('theta: %5.2f', control.theta);
            strLambda = sprintf('lambda: %5.2f', control.lambda);
            strBundle = sprintf('bundle size: %d', control.bundle_limit);
            strLb_mode = sprintf('LB mode [0/1]: %d', control.lb_mode);
%            strSubSolver = sprintf('Subproblem solver [0/1]: %d', control.SubSolver);
%            strSubITER  = sprintf('Subproblem iteration: %d', control.SubITER);
%            strSubEPS = sprintf('Subproblem epsilon: %.2e', control.SubEPS);
%            strSubLambda = sprintf('Subproblem lambda: %5.2f', control.SubLambda);
             
            irun=menu(strP, strIter, strEpsilon, strTheta, strLambda, strBundle, strLb_mode, ...
            'Run','->');
            if irun == 8,
                break;
            end;
            if irun == 1,
                control.iter_limit = input('iteration limit > ');
            end;
            if irun == 2,
                control.epsilon = input('epsilon > ');
            end;
            if irun == 3,
                control.theta = input('theta > ');
            end;
            if irun == 4,
                control.lambda = input('lambda > ');
            end;
            if irun == 5,
                control.bundle_limit = input('bundle size > ');
            end;
            if irun == 6,
                control.lb_mode = input('LB mode > ');
            end;          
            if irun == 7,
                frep = fopen('report.txt','a');
                fprintf(frep,'** run APL for: %s**\n',name);
                fclose(frep);
                
                flag=input(sprintf('To run APL_SP_uniform[%1d,%1d] [y/n] > ',control.iter_limit,control.bundle_limit),'s');
                if flag(1)=='y',
                    APL_SP_uniform(control, domain, data);
                end;
                
                flag=input(sprintf('To run NERML_SP[%1d,%1d] [y/n] > ',control.iter_limit,control.bundle_limit),'s');
                if flag(1)=='y',
                    NERML_SP(control, domain, data);
                end;
                
                flag=input(sprintf('To run APL_SP_Nonsmooth[%1d,%1d] [y/n] > ',control.iter_limit,control.bundle_limit),'s');
                if flag(1)=='y',
                    APL_SP_nonsmooth(control, domain, data);
                end
            end;
        end;
    end;
end;





