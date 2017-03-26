function [sol] = ProxMappingCVX(domain, Bundle, data, BarX, f, g, x_lbt, c, LS)

x = zeros(domain.n,1);
delta = 1.0e-16 / domain.n;
domega_c = 1 + log(c + delta);

if Bundle.size == 0,
    cvx_begin
    cvx_quiet(true);
        variable x(domain.n);
        minimize -sum(entr(x+delta) -domega_c.*x);
        subject to
            x == simplex(domain.n);
            BarX.a' * x >= BarX.b;
            f + g' * (x - x_lbt) <= LS;
    cvx_end
else,
    cvx_begin
    cvx_quiet(true);
        variable x(domain.n);
        minimize -sum(entr(x+delta) -domega_c.*x);
        subject to
            x == simplex(domain.n);
            BarX.a' * x >= BarX.b;
%              f + g' * (x - x_lbt) <= LS;
            Bundle.const(1:Bundle.size) + Bundle.matrix(1:Bundle.size, :) * x <= LS * ones(Bundle.size, 1);
    cvx_end
end;


status = sprintf(cvx_status);
if strcmp(status,'Solved') ~= 1
    disp(status);    
end
sol = x;