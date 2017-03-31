function [loss , obj , grad_1, grad_2] = evaluation_oracle2(zz,data)

    %% Non-convex S3VM problem
    %     er1 = data.test_matrix*zz;
    %     obj1 = tanh(data.test_lable.*er1);
    %     loss1 =  sum(ones(data.vali,1) - obj1)/data.vali;
    %     obj = loss1 + data.lambda*norm(zz)^2;
    %     loss = 100*sum(abs(sign(er1)-data.test_lable))/(2*data.vali);
    %     grad = data.test_matrix'*(data.test_lable.*(obj1.^2-ones(data.vali,1)))/data.vali+ 2*data.lambda*zz;
    obj1 = data.lambda2*mean(max(0,1-data.test_lable.*(data.test_matrix(1:data.vali,:)*zz)).^2);
    obj2 = data.lambda3*mean(exp(-5*(data.test_matrix(data.vali+1:end,:)*zz).^2));
    obj = 0.5*data.lambda1*norm(zz)^2 + obj1+obj2;
    loss =obj;
    GR_b = 1-data.test_lable.*(data.test_matrix(1:data.vali,:)*zz);
    GR_b = 2*GR_b.*(GR_b>0);
    GR2 = [mean(-data.test_lable.*GR_b);data.test_matrix(1:data.vali,2:end)'*(-data.test_lable.*GR_b)/data.vali];
    GR_b = -10*(data.test_matrix(data.vali+1:end,:)*zz).*exp(-5*(data.test_matrix(data.vali+1:end,:)*zz).^2);
    GR3 = [mean(GR_b);data.test_matrix(data.vali+1:end,2:end)'*(GR_b)/data.vali];
    grad_1 = data.lambda1*zz+data.lambda2*GR2+ data.lambda3*GR3;
    grad_2 = grad_1;
    a = zz(1)-data.gamma*grad_1(1);
    if a > 2*data.ratio-1+data.toler
        a = 2*data.ratio-1+data.toler;
    elseif a < 2*data.ratio-1-data.toler
        a = 2*data.ratio-1-data.toler;
    end  
    grad_2(1) = (zz(1)-a)/data.gamma;
    
