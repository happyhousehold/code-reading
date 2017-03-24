%function SG = SGradient(x,y,z_md,lamda)
GR = -2*(y(index) - dot(z_md, x(index,:)))*x(index,:)'+ 2*lamda*z_md;
