% Generating or Loading initila sample
while  running == 0
    data_name = ['Initial Sample : ', filename];
    mtitle = data_name;
    awr = menu(mtitle, 'Generating new sample' , 'Load sample ' , 'Run Algorithms' , ' Quit');
    %awr = menu(mtitle, 'Load training sample ' , 'load testing sample' , 'Run Algorithms' , ' Quit');
    if awr == 1
        data.dim = input('Enter the dimension of problem : ');
        data.size = input('Enter the size of data set : ');
        data.spr = input('Enter the sparsity of data set : '); 
        data.seed = input('Enter the seed for generating random numbers : ');
        instance_matrix = ceil(sprand(data.size,data.dim, data.spr));
        z_sample = 2*rand(data.dim, 1)-ones(data.dim,1);
        data.sample = z_sample;
        label_vector = sign(instance_matrix*z_sample);
        d = data.dim;
        sz = data.size;
        estimate_matrix = instance_matrix(1:N_initial,:);
        estimate_lable = label_vector(1:N_initial);
        data.test_matrix = instance_matrix(N_training+1:end ,:);
        data.test_lable = label_vector(N_training+1: end,:);
        z_ini = 2*rand(data.dim, 1)-ones(data.dim,1);
        data.vali = length(data.test_lable);
        filename = ['SVM','-n',num2str(d)];
        writing_seed;
        save (filename, 'data', 'z_ini' , 'z_sample', 'estimate_matrix', 'estimate_lable');
        clear instance_matrix label_vector
    elseif awr == 2
        [filename,PathName,FilterIndex] = uigetfile('*.mat');
        [test_lable, test_matrix] = libsvmread(filename);
        test_matrix = sparse(test_matrix);
    elseif  awr == 3
        if filename == ' '
            error('You should choose an initial sample');
        end
        running = 1;
    elseif  awr == 4
        break;
    end
end
if  awr == 4
    break;
end



