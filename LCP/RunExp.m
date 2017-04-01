%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MATLAB Code for the Projection-free Gradient (PFG) Methods            %
%     Author: Guanghui (George) Lan                                      %
%     Institute: University of Florida, Industrial & Systems Engineering %
%     @All rights reserved 2013                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this file implements an interface for testing the PFG algorithms %
%clear all;
name='noname';
fname='noname.ins';

flag_dat = 0;

% control parameters
control.theta = 0.5;
control.lambda = 0.5;
control.iter_limit = 1000;
control.epsilon = 1.0e-6;
control.bundle_limit = 30;
control.lb_mode = 0; %% 0 if minimize a linear function, 1 to minimize the
                     %% piecewise linear function.
                     
while(1==1)
    strD=sprintf('Data: %s',name);
    imenu=menu(strD,'Generate data','Load data', 'Run', 'View protocol','Quit');
    if imenu==4,
        edit report.txt
    end;
    if imenu==5,
        break;
    end;
    if imenu == 1,
        % default settings
        % Generate quadratic programming instances.
        dataF.domain.type = 1;
        dataF.domain.n = 2000;
        dataF.domain.m = 100;
        dataF.domain.R = 10;
        dataF.data.mode = 2;
   %     dataF.data.cn   = 1e3;
        dataF.data.datasparsity = 0.5;
        dataF.data.solsparsity = 0.2;
        dataF.data.rank = 0;
        dataF.data.noice = 0;

        dataF.domain.type = input('domain type [0/1/2/3/4 for hyb/box/simplex/ball/Spectahedron] > ');
        dataF.domain.n = input('number of variables > ');
        dataF.domain.m = input('number of rows > ');

        dataF.domain.R = 1; 
        if dataF.domain.type == 0
           dataF.domain.R = input('radious > ');
        end
        %if dataF.domain.type == 4
        %    dataF.data.rank = input('data rank > ');
        %end
  %      dataF.domain.cn  = input('condition number > ');
        dataF.data.datasparsity  = input('data density > ');
   %     dataF.data.solsparsity  = input('solution density > ');
    %    dataF.data.noice = input('noice level > ');
        
        dataF.data = NewInst(dataF.domain.type, dataF.domain.m, dataF.domain.n,...
            dataF.domain.R, dataF.data.datasparsity, dataF.data.solsparsity, dataF.data.noice);
                
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
        control.iter_limit = 300;
        control.optstepsize = 0;
        %%%%% Run 
        while(1==1)
            
            % choose options
            strP = sprintf('Processing data %s',name);
            strIter = sprintf('Iteration limit: %d', control.iter_limit);
            strOptStep = sprintf('Optimal stepsize [0/1]: %d', control.optstepsize);
             
            irun=menu(strP, strIter, strOptStep, 'Run','->');
            if irun == 4,
                break;
            end;
            if irun == 1,
                control.iter_limit = input('iteration limit > ');
            end;
            if irun == 2,
                control.optstepsize = input('use optimal stepsize > ');
            end;           
          
            if irun == 3,
                frep = fopen('report.txt','a');
                fprintf(frep,'** run projection-free methods for: %s**\n',name);
                fclose(frep);
                               
                %%% compute a random starting point.
                x0 = GetRandInitPt(dataF.domain.type, dataF.domain.n, dataF.domain.R);
                [fx0, gx0] = FirstOrderOracleQP(dataF.data,dataF.domain,x0);
                frep = fopen('report.txt','a');
                str = sprintf('Initial value: %.6e\n', fx0);
                disp(str);
                fprintf(frep,str);
                fclose(frep);

                
                flag=input(sprintf('To run CG[%1d,%1d] [y/n] > ',control.iter_limit,control.optstepsize),'s');
                if flag(1)=='y',
                    frep = fopen('report.txt','a');
                    fprintf(frep,sprintf('To run CG[%1d,%1d]',control.iter_limit,control.optstepsize));
                    fclose(frep);
                    CG(control, dataF.domain, dataF.data, x0);
                end;
                
                flag=input(sprintf('To run PA_CG[%1d,%1d] [y/n] > ',control.iter_limit,control.optstepsize),'s');
                if flag(1)=='y',
                    frep = fopen('report.txt','a');
                    fprintf(frep,sprintf('To run PA_CG[%1d,%1d]',control.iter_limit,control.optstepsize));
                    fclose(frep);
                    PA_CG(control, dataF.domain, dataF.data, x0);
                end;
                
                flag=input(sprintf('To run PDA_CG[%1d,%1d] [y/n] > ',control.iter_limit,control.optstepsize),'s');
                if flag(1)=='y',
                    frep = fopen('report.txt','a');
                    fprintf(frep,sprintf('To run PDA_CG[%1d,%1d]',control.iter_limit,control.optstepsize));
                    fclose(frep);
                    PDA_CG(control, dataF.domain, dataF.data, x0);
                end;
                
    %            flag=input(sprintf('To run PDA_Restart_CG[%1d,%1d] [y/n] > ',control.iter_limit,control.optstepsize),'s');
    %            if flag(1)=='y',
    %                PDA_Restart_CG(control, dataF.domain, dataF.data, x0);
    %            end;
                
                if dataF.domain.type == 1 %% box constraint
                    flag=input(sprintf('To run NEST_QP_Original[%1d] [y/n] > ',control.iter_limit),'s');
                    if flag(1)=='y',
                        NEST_QP_Original(control, dataF.domain, dataF.data, x0);
                    end
                end
                
                if dataF.domain.type == 2 %% simplex constraint
                    
                    flag=input(sprintf('To run NEST_QP_Original_Simplex_Eu[%1d] [y/n] > ', control.iter_limit),'s');
                    if flag(1)=='y',
                        NEST_QP_Original_Simplex_Eu(control, dataF.domain, dataF.data, x0);
                    end
                end
                
                if dataF.domain.type == 4 %% Spectahedron constraint
                    flag=input(sprintf('To run NEST_QP_Original_Spect_Eu[%1d] [y/n] > ',control.iter_limit),'s');
                    if flag(1)=='y',                       
                        NEST_QP_Original_Spect_Eu(control, dataF.domain, dataF.data, x0);
                    end
                end
            end;
        end;
   end;  
end;





