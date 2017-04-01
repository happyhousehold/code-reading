% compute lambdmax(A A')
function [v,lmax] = myEig(A, m, tol)

v = rand(m,1);

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
