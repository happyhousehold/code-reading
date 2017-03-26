function [lb_t] = ComputeNewLBCVX(domain, Bundle, data, BarX, f, g, x_lb, LS, mode)
x = zeros(domain.n,1);
t = 0;

if Bundle.size == 0,
    cvx_begin
    cvx_quiet(true);
        variable x(domain.n);
        minimize (f + g' * (x - x_lb));
        subject to
            x == simplex(domain.n);
            BarX.a' * x >= BarX.b;
    cvx_end
else,
    if mode == 0,
        cvx_begin
        cvx_quiet(true);
            variable x(domain.n);
            minimize (f + g' * (x - x_lb));
            subject to
                x == simplex(domain.n);
                BarX.a' * x >= BarX.b;
                Bundle.const(1:Bundle.size) + Bundle.matrix(1:Bundle.size,:) * x <= LS * ones(Bundle.size, 1);
        cvx_end

    else,
        cvx_begin
        cvx_quiet(true);
            variables x(domain.n) t;
            minimize t;
            subject to
                x == simplex(domain.n);
                BarX.a' * x >= BarX.b;
                Bundle.const(1:Bundle.size) + Bundle.matrix(1:Bundle.size,:) * x <= LS * ones(Bundle.size, 1);
                f + g' * (x - x_lb) < t;
                Bundle.const(1:Bundle.size) + Bundle.matrix(1:Bundle.size,:) * x < t;
        cvx_end
    end
end;

status = sprintf(cvx_status);
if strcmp(status,'Solved') ~= 1
    disp(status);    
end
lb_t = cvx_optval;