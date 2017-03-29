% Generating or Loading initila sample
while  running == 0
    data_name = ['Instance name : ', filename];
    mtitle = data_name;
    awr = menu(mtitle, 'Load sample ' , 'Generating new sample' , 'Run Algorithms' , ' Quit');
    if awr == 1
        %[filename,PathName,FilterIndex] = uigetfile('*.mat');
        if  strcmp(data.problem, 'Least_sqr') == 1
            dir *Least_sqr*.mat
        elseif  strcmp(data.problem, 'Reg') == 1
            dir *Reg*.mat
        elseif  strcmp(data.problem, 'SVM') == 1
            dir *SVM*.mat
        elseif   strcmp(data.problem, 'Trig') == 1
            dir *Trig*.mat
        end
        filename = input('Enter the name of instance: ');
        load (filename);
    elseif awr == 2
        fprintf('  \n');
        d = input('Enter the dimension of problem : ');
        data.dim = d;
        data.st = input('Enter the standard variation of noise : ');
        st = data.st;
        data.spr = input('Enter the sparsity of data : ');
        data.seed = input('Enter the seed for generating random numbers : ');
        er_initial = normrnd(0,data.st,N_initial,1);
        x_initial = sprand(N_initial, data.dim, data.spr);
        z_sample = rand(data.dim, 1);
        z_ini = rand(data.dim, 1);
        y_initial = x_initial*z_sample + er_initial;
        sparsity = 100*data.spr;
        data.validx = sprand(data.vali, data.dim, data.spr);
        data.validy =  data.validx*z_sample + normrnd(0,data.st,data.vali,1);
        std = 10*st;
        filename = ['Reg-n',num2str(d),'-st',num2str(std)];
        writing_seed;
        save (filename, 'N_initial','data', 'er_initial', 'x_initial', 'z_sample', 'z_ini', 'y_initial');
    elseif  awr == 3
        if filename == ' '
            %error('You should choose an initial sample');
            msgbox('Please choose an initial sample before running algorithms!','Error','error');
        else
            running = 1;
        end
    elseif  awr == 4
        break;
    end
end
if  awr == 4
    break;
end



