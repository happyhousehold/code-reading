% Generating or Loading initila sample
while  running == 0
    data_name = ['Initial Sample : ', filename];
    mtitle = data_name;
    awr = menu(mtitle, 'Load sample ' , 'Generating new sample' , 'Run Algorithms' , ' Quit');
    if awr == 1
        [filename,PathName,FilterIndex] = uigetfile('*.mat');
        load (filename);
    elseif awr == 2
        fprintf('  \n');
        N_initial = input('Enter Number of initial samples : ');
        d = input('Enter the dimension of problem : ');
        st = input('Enter the standard variation of noise : ');
        er_initial = normrnd(0,st,1,N_initial)';
        x_initial = rand(N_initial,d);
        % z is the substitution for x in the algorithm
        z_initial = rand(1,d)';
        z_ini = rand(1,d)';
        y_initial = x_initial*z_initial + er_initial;
        M = 0;
        epsln = 0.01;
        filename = ['N', num2str(N_initial),'-d',num2str(d),'-st', num2str(st)];
        save (filename, 'N_initial','d', 'st', 'er_initial', 'x_initial', 'z_initial', 'z_ini', 'y_initial', 'M' ,'epsln'); 
    elseif  awr == 3
        if filename == ' '
            error('You should choose an initial sample');
        end
        running = 1;
    elseif  awr == 4
        break;
    end
end



