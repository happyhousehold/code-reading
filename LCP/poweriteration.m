% compute lambdamax(A A')
function lmax = poweriteration(A, m, n, tol)

v = rand(m,1);

v = v / norm(v);

while 1==1
    vbar = A * A' * v;
    vbar = vbar / norm(vbar);
    
    acc = vbar' * v;
    
    if acc > 1 - tol,
        break;
    end
    
    v = vbar;
    
end

lmax = v' * A* A' * v;