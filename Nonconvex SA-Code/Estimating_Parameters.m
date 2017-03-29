%% Estimating problem parameters
initial_solution = R_ini * z_ini;

if (strcmp(data.problem, 'Least_sqr') == 1) || (strcmp(data.problem, 'Reg') == 1)
    %% Non-convex Regression
    sig = zeros(1, 10);
    for i = 1 : 10
        initial = rand(data.dim, 1);
        stch_grad = zeros(20,data.dim);
        for j = 1 : 20
            stch_grad(j,:)= (2*(dot(initial, x_initial((i-1)*20+j,:))-y_initial((i-1)*20+j))*x_initial((i-1)*20+j,:)'+ data.lambda*p*initial.*(initial.^2+ep^2*ones(data.dim,1)).^(p/2-1))';
        end
        sig(i) = sqrt(sum(var(stch_grad)));
    end
    sigma = max(sig);
    A = 2*(x_initial'*x_initial)/N_initial;
    L = real(max(eig(A)))+ data.lambda*p*ep^(p-2);
    data.sample = z_sample;
    data.eps = ep;
    data.p = p;
    d = data.dim;
    [loss_0 , e_0 , g_0] = evaluation_oracle(initial_solution,data);
    
elseif (strcmp(data.problem, 'SVM') == 1)
    %% Non-convex SVM problem
    sig = zeros(1,10);
    A = zeros(data.dim);
    for i =1 : N_initial
        A = A + estimate_matrix(i,:)'*estimate_matrix(i,:);
    end
    A = 2*A/N_initial+2*data.lambda*eye(data.dim);
    L= max(real(eig(A)));
    sig = zeros(1, 10);
    for i = 1 : 10
        initial = rand(data.dim, 1) - rand(data.dim, 1);
        stch_grad = zeros(20,data.dim);
        for j = 1 : 20
            stch_grad(j,:) = estimate_lable((i-1)*20+j)*(tanh(estimate_lable((i-1)*20+j)*dot(initial, estimate_matrix((i-1)*20+j,:)'))^2-1)*estimate_matrix((i-1)*20+j,:)' + 2*data.lambda*initial;
        end
        sig(i) = sqrt(sum(var(stch_grad)));
    end
    sigma = max(sig);
    [loss_0, e_0 , g_0] = evaluation_oracle(initial_solution,data);
    
elseif (strcmp(data.problem, 'Trig') == 1)
    %% Trigonometric functions
    sigma = sqrt(data.dim)*data.st;
    for i = 1 : 10
        x1 = rand(data.dim,1);
        x2 = rand(data.dim,1);
        G1 = zeros(data.dim,1);
        G2 = zeros(data.dim,1);
        L1 = zeros(10,1);
        for j=1:20
            z = x1;
            SGradient;
            G1 = G1+ GR;
            z = x2;
            SGradient;
            G2 = G2+GR;
        end
        G1 = G1/20;
        G2 = G2/20;
        L1(i) = norm(G1 - G2)/norm(x1 - x2);
    end
    L = max(L1);
    [loss_0 , e_0 , g_0] = evaluation_oracle(initial_solution,data);
    
end