
function [loss , obj , grad_1, grad_2] = evaluation_oracle(zz,data)

%% Non-convex S3VM problem
data.test_matrix2 = [ones(2*data.vali,1),sprandn(2*data.vali,data.dim-1,data.spr)];
data.test_lable2 = sign(data.test_matrix2(1:data.vali,:)*data.sample);
obj1 = data.lambda2*mean(max(0,1-data.test_lable2.*(data.test_matrix2(1:data.vali,:)*zz)).^2);
obj2 = data.lambda3*mean(exp(-5*(data.test_matrix2(data.vali+1:end,:)*zz).^2));
obj = 0.5*data.lambda1*norm(zz)^2 + obj1+obj2;
loss =obj;
%grad = data.test_matrix2'*(data.test_lable2.*(obj1.^2-ones(data.vali,1)))/data.vali+ 2*data.lambda*zz;
%Lasso
%     GG=data.test_matrix2'*(data.test_lable2.*(obj1.^2-ones(data.vali,1)))/data.vali;
%     grad = (zz- sign(zz - data.gamma*GG).*max(abs(zz - data.gamma*GG)-data.gamma*data.lambda,0))/data.gamma;
%     grad2  = (abs(zz)>0.0001).*(GG+data.lambda*sign(zz));
GR_b = 1-data.test_lable2 .*(data.test_matrix2(1:data.vali,:)*zz);
GR_b = 2*GR_b.*(GR_b>0);
GR2 = [mean(-data.test_lable2.*GR_b);data.test_matrix2(1:data.vali,2:end)'*(-data.test_lable2.*GR_b)/data.vali];
GR_b = -10*(data.test_matrix2(data.vali+1:end,:)*zz).*exp(-5*(data.test_matrix2(data.vali+1:end,:)*zz).^2);
GR3 = [mean(GR_b);data.test_matrix2(data.vali+1:end,2:end)'*(GR_b)/data.vali];
grad_1 = data.lambda1*zz+data.lambda2*GR2+ data.lambda3*GR3;
grad_2 = grad_1;

