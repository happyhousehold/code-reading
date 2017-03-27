% compute lambdamax(A)
function [lmax,v] = poweriteration(A, n, tol)

v = rand(n,1);

v = v / norm(v);

while 1==1
    vbar = A * v;
    vbar = vbar / norm(vbar);
    
    acc = vbar' * v;
    
    if acc > 1 - tol,
        break;
    end
    
    v = vbar;
    
end

lmax = v' * A * v;