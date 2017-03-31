% Generating or Loading initila sample
while  running == 0
    data_name = ['Instance name : ', filename];
    mtitle = data_name;
    awr = menu(mtitle, 'Load sample ' , 'Generating new sample' , 'Run Algorithms' , ' Quit');
    if awr == 1
        dir *S3VM*.mat
        filename = input('Enter the name of instance: ' , 's');
        load (filename);
    elseif awr == 2
        load last_seed.txt
        last_seed = last_seed +10000;
        delete last_seed.txt
        save -ascii last_seed.txt last_seed
        Sample_generation_S3VM; 
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
