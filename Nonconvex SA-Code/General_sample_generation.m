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
        filename = input('Enter the name of instance: ' , 's');
        load (filename);
    elseif awr == 2
        load last_seed.txt
        last_seed = last_seed +10000;
        delete last_seed.txt
        save -ascii last_seed.txt last_seed
        if  strcmp(data.problem, 'Least_sqr') == 1
            Sample_generation_Reg;
        elseif  strcmp(data.problem, 'Reg') == 1
            Sample_generation_Reg;
        elseif  strcmp(data.problem, 'SVM') == 1
            Sample_generation_SVM;
        elseif   strcmp(data.problem, 'Trig') == 1
            Sample_generation_Trig;
        end
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
