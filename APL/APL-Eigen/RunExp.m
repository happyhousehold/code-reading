%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MATLAB Code for the Accelerated Prox-Level (APL) algorithm            %
%     Author: Guanghui (George) Lan                                      %
%     Institute: University of Florida, Industrial & Systems Engineering %
%     @All rights reserved 2010                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this file implements an interface for testing the APL algorithm %
%clear all;
name='noname';
fname='noname.ins';

flag_dat = 0;

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
    strD=sprintf('Data: %s',name);
    imenu=menu(strD,'Generate data','Load data','Run APL', ...
        'Run NEST', 'Run IPM', 'View protocol','Quit');
    if imenu==6,
        edit report.txt
    end;
    if imenu==7,
        break;
    end;
    if imenu == 1,
        % default settings    
        dataF.domain.n = 100;
        dataF.domain.m = 50;
        dataF.domain.type = 'se';
            
        dataF.domain.n = input('number of variables > ');
        dataF.domain.m = input('dimension of matrices > ');         
        dataF.data.density  = input('density > ');
        dataF.data.diagonal = input('diagonal > ');

        dataF.data = NewEigenInst(dataF.domain.n, dataF.domain.m, dataF.data.density, dataF.data.diagonal);
        dataF.domain.center = zeros(dataF.domain.n,1);
        dataF.domain.axes = ones(dataF.domain.n,1);

        name = input('dataname > ','s');
        frep = fopen('report.txt','a');
        fprintf(frep,'** generate data: %s**\n',name);
        fclose(frep);
        str = sprintf('%s=dataF;',name);
        eval(str);
        str = sprintf('save %s.mat %s',name,name);
        eval(str);
        flag_dat = 1;
    end;
    
    if imenu == 2,
        ls *.mat
        name = input('name > ','s');
        frep = fopen('report.txt','a');
        fprintf(frep,'** load data: %s\n',name);
        fclose(frep);
        str = sprintf('load %s.mat',name);
        eval(str);
        str = sprintf('dataF=%s;',name);
        eval(str);
        flag_dat = 1;
    end;
    
    %%%%%%%%%
    if imenu==3,
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
        %% the following options are valid only if control.SubSolver = 1.
        control.SubITER = 200;
        control.SubEPS = 1e-10;
        control.SubLambda = 0.5;
        
       
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
            strDType = sprintf('DomainType [0(simplex)/>=1(box)]:%d', control.DomainType);
%            strSubSolver = sprintf('Subproblem solver [0/1]: %d', control.SubSolver);
%            strSubITER  = sprintf('Subproblem iteration: %d', control.SubITER);
%            strSubEPS = sprintf('Subproblem epsilon: %.2e', control.SubEPS);
%            strSubLambda = sprintf('Subproblem lambda: %5.2f', control.SubLambda);
             
            irun=menu(strP, strIter, strEpsilon, strTheta, strLambda, strBundle, strLb_mode, strDType, ...
            'Run','->');
            if irun == 9,
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
                control.DomainType = input('Domain type > ');
            end;
            if irun == 8,
                frep = fopen('report.txt','a');
                fprintf(frep,'** run APL for: %s**\n',name);
                fclose(frep);
                dataF.data.DomainType = control.DomainType;
                
                flag=input(sprintf('To run APL_Eigen_basic[%1d,%1d] [y/n] > ',control.iter_limit,control.bundle_limit),'s');
                if flag(1)=='y',
                    APL_Eigen_basic(control, dataF.domain, dataF.data);
                end;
                
      %          flag=input(sprintf('To run APL_Eigen_basic_fixsmoothing[%1d,%1d] [y/n] > ',control.iter_limit,control.bundle_limit),'s');
      %          if flag(1)=='y',
      %              APL_Eigen_basic_fixsmoothing(control, dataF.domain, dataF.data);
      %          end;
                
                flag=input(sprintf('To run APL_Eigen_uniform[%1d,%1d] [y/n] > ',control.iter_limit,control.bundle_limit),'s');
                if flag(1)=='y',
                    APL_Eigen_uniform(control, dataF.domain, dataF.data);
                end;
                
                flag=input(sprintf('To run APL_Eigen_basic_nonsmooth[%1d,%1d] [y/n] > ',control.iter_limit,control.bundle_limit),'s');
                if flag(1)=='y',
                    APL_Eigen_basic_nonsmooth(control, dataF.domain, dataF.data);
                end;
                
                flag=input(sprintf('To run NERML_Eigen[%1d,%1d] [y/n] > ',control.iter_limit,control.bundle_limit),'s');
                if flag(1)=='y',
                    NERML_Eigen(control, dataF.domain, dataF.data);
                end;
            end;
        end;
    end;
   
   if imenu==4,
        if flag_dat==0,
            error('No data specified');
            break;
        end;

        % default control options                         
        control.iter_limit = 1000;
 
        %%%%% Run 
        control.stepfactor = 1.0;
        while(1==1)
            
            % choose options
            
            strP = sprintf('Processing data %s',name);
            strIter = sprintf('iteration limit: %d', control.iter_limit);
            strEpsilon = sprintf('epsilon: %.2e', control.epsilon);
            strDType = sprintf('DomainType [0(simplex)/1(box)]:%d', control.DomainType);
            strFactor = sprintf('step factor: %.2e', control.stepfactor);

            irun=menu(strP, strIter, strEpsilon, strDType, strFactor, 'Run','->');
            if irun == 6,
                break;
            end;
            if irun == 1,
                control.iter_limit = input('iteration limit > ');
            end;
            
            if irun == 2,
                control.epsilon = input('epsilon > ');
            end;
            
            if irun == 3,
                control.DomainType = input('DomainType [0(simplex)/1(box)] > ');
            end;
            
            if irun == 4,
                control.stepfactor = input('step factor > ');
            end;
            
            if irun == 5,  
                frep = fopen('report.txt','a');
                fprintf(frep,'** run NEST for: %s**\n',name);
                fclose(frep);
                start_NEST=clock;
                dataF.data.DomainType = control.DomainType;
%                NEST_QP(control, dataF.domain, dataF.data);
                if dataF.data.DomainType == 0,
                    dataF.domain.type = 'se';
                    NEST_Eigen(control, dataF.domain, dataF.data);
                else,
                    flag=input(sprintf('To run NEST [y/n] > '),'s');
                    if flag(1)=='y',
                        NEST_Eigen_box(control, dataF.domain, dataF.data);
                    end;
                    
                    flag=input(sprintf('To run NEST_Search [y/n] > '),'s');
                    if flag(1)=='y',
                        NEST_Eigen_box_Search(control, dataF.domain, dataF.data);
                    end;
                end;
                time_NEST=etime(clock,start_NEST);
                disp(sprintf('time=%5.2f', time_NEST));
            end;
        end;
    end;
    
    if imenu==5,
        if flag_dat==0,
            error('No data specified');
            break;
        end;
        
        start_CVX=clock;
        cvx_begin
            cvx_quiet(false);
                variable x(dataF.domain.n+1);
                minimize lambda_max(reshape(dataF.data.exA * x, dataF.data.m, dataF.data.m))
                subject to
                    x(1) == 1;
                    x(2:dataF.domain.n+1) == simplex(dataF.domain.n) ;                  
        cvx_end
        
        time_CVX=etime(clock,start_CVX);
        
        
        disp(sprintf('UB=%.6e, time=%5.2f', cvx_optval, time_CVX));
        xopt = x(2:dataF.domain.n+1);
        str = sprintf('save %s_sol.mat xopt',name);
        eval(str);
    end;
end;





