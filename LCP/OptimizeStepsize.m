function [astar] = OptimizeStepsize(data,domain,xt, yt);

%%%% bisection

l = 0;
u = 1;

delta = xt - yt;
while (u - l > 1e-6)
    testa = (l + u)/2;
    [f,g] = FirstOrderOracleQP(data,domain, yt + testa * delta);
    dir = g' * delta;
    if dir > 0,
        u = testa;
    else
        l = testa;
    end
end

astar = u;