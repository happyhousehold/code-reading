function ini = GetRandInitPt(type, n, R)

%%%% determine the number of nonzeros in the solution
%%%% dependent on the domain type, generate an optimal solution

if type == 0 %box intersected with simplex
    ini = [rand(n,1)];
    l1 = sum(ini);
    
    if l1 > R
        ini = ini / l1 * R;
    end
end

if type == 1 %% box
    ini = [R*(rand(n,1))];   
end

if type == 2  %%%% simplex 
%    x=[sort(rand(n-1,1));1];
%    y=[0;x(1:n-1)];
%    ini=R * (x-y);
    ini = R * ones(n,1)/n;
end

if type == 3%%% ball
    tx = [rand(n,1)];
    ini = R * tx / norm(tx);
end

if type == 4 %%% spectahedron
%    MF = sprand(n, n, 0.5);
%    sol = sparse(MF * MF');
%    tr = trace(sol);
%    ini = reshape(sol/tr*R,n^2,1);    
    ini = reshape(eye(n)*R/n, n^2,1);
end



