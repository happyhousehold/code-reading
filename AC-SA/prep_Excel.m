objective_final = mean(objective_value);
var_objective_final = var (objective_value);
run_time = mean(r_time);
mean_data = {'runtime' ,'objfun_final' , ' variance'; run_time , objective_final, var_objective_final}
mean_data2 = ['runtime' ,'objfun_final' , ' variance'];
sheet_name = ['lamda-', num2str(lamda)];
%data_showing;