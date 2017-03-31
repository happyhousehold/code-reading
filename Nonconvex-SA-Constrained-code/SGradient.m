%% Non-convex S3VM PROBLEM
Data_generation;
%GR = y*(tanh(y*dot(z,x))^2-1)*x+ 2*data.lambda*z;
%GR = x'*(y.*(tanh(y.*(x*z)).^2-1))/data.m;
%GR_b = y<x,w>+b
GR_b = 1-y.*(x1*z);
GR_b = 2*GR_b.*(GR_b>0);
GR2 = [mean(-y.*GR_b);x1(:,2:end)'*(-y.*GR_b)/m];
GR_b = -10*(x2*z).*exp(-5*(x2*z).^2);
GR3 = [mean(GR_b);x2(:,2:end)'*(GR_b)/m];
GR =data.lambda1*z+data.lambda2*GR2+ data.lambda3*GR3;
