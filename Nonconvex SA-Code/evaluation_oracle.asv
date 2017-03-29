function [loss , obj , grad] = evaluation_oracle(zz,data)

if (strcmp(data.problem, 'Least_sqr') == 1) || (strcmp(data.problem, 'Reg') == 1)
    %% Non-convex Regression
    data.valx = sprand(data.vali, data.dim, data.spr);
    data.valy = data.valx*data.sample + normrnd(0,data.st,data.vali,1);
    loss = (1/data.vali)*norm(data.valx*zz - data.valy)^2;
    obj = loss + data.lambda*sum((zz.^2+data.eps^2*ones(data.dim,1)).^(data.p/2));
    loss = obj;
    grad = (2/data.vali)*data.valx'*(data.valx*zz - data.valy)+ data.lambda*data.p*zz.*(zz.^2+data.eps^2*ones(data.dim,1)).^(data.p/2 -1);
    
elseif (strcmp(data.problem, 'SVM') == 1)
    %% Non-convex SVM problem
    data.test_matrix2 = ceil(sprand(data.vali,data.dim, data.spr));
    data.test_lable2 = sign(data.test_matrix2*data.sample);
    er1 = data.test_matrix2*zz;
    obj1 = tanh(data.test_lable2.*er1);
    loss1 =  sum(ones(data.vali,1) - obj1)/data.vali;
    obj = loss1 + data.lambda*norm(zz)^2;
    loss = 100*sum(abs(sign(er1)-data.test_lable2))/(2*data.vali);
    grad = data.test_matrix2'*(data.test_lable2.*(obj1.^2-ones(data.vali,1)))/data.vali+ 2*data.lambda*zz;
    
elseif (strcmp(data.problem, 'Trig') == 1)
    %% trigonometric functions
    err = normrnd(0, data.st, data.vali,data.dim);
    obj1 = (data.dim - sum(cos(zz)))*ones(data.dim, 1)+cumsum(ones(data.dim, 1)).*(ones(data.dim, 1)-cos(zz))-sin(zz);
    obj = norm(obj1)^2;
    loss = obj;
    grad1 = sum(obj1)*sin(zz)+obj1.*(cumsum(ones(data.dim, 1)).*sin(zz)-cos(zz));
    grad = 2*grad1 + mean(err)';
end