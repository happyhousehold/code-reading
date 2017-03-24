objective_final = mean(objective_value);
run_time = mean(r_time);
final_data = [run_time , objective_final];
fprintf ('Runtime:% 4.2f     Final objective value:% 6.2f\n', run_time , objective_final);

