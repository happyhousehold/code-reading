gradient_final = norm(gg)^2;
final_data = [objective_value, final_loss];
%fprintf ('Runtime:% 4.2f     Final loss:% 7.4f\n  Final objective value:% 7.4f\n    Final gradient:% 7.4f\n', run_time , final_loss, objective_value, gradient_final);
