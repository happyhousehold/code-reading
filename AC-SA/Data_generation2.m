x = rand(N_iter,d);
er = normrnd(0,st,1,N_iter)';
y = x*z_sample + er;