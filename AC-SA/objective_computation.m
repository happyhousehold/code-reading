clear x y er
r_time(j) = toc;
x_vali = rand(N_vali,d);
y_vali = x_vali*z_sample + normrnd(0,st,1,N_vali)';
objective_value(j) = (1/N_vali)*(y_vali - x_vali*final_solution')'*(y_vali - x_vali*final_solution')+ lamda*norm(final_solution)^2;

