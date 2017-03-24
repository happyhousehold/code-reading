% data1 = [mean_data{1}, '     ', mean_data{3}, '     ', mean_data{5}, '   ', mean_data{7}];
% blnk = '                                                                            ';   
% data2 = [num2str(mean_data{2}), '      ' , num2str(mean_data{4}), '              ' ,num2str(mean_data{6}), '              ' ,num2str(mean_data{8})];
% data_show = {blnk; data1 ; blnk; blnk; data2; blnk; blnk; blnk};
data1 = [mean_data{1}, '     ', mean_data{3},'     ', mean_data{5}];
blnk = '                                                          ';   
data2 = [num2str(mean_data{2}), '      ' , num2str(mean_data{4}),'      ' , num2str(mean_data{6})];
data_show = {blnk; data1 ; blnk; blnk; data2; blnk; blnk};
msgbox(data_show);
