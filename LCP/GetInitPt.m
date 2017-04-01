function ini = GetInitPt(type, n, R)

%%%% determine the number of nonzeros in the solution
%%%% dependent on the domain type, generate an optimal solution
if type == 1 %% box
    ini = zeros(n,1);    
else
    if type == 2  %%%% simplex 
        
        ini=R * [1;zeros(n-1,1)]
    else 
        if type == 3 %% ball
            ini = R * [1;zeros(n-1,1)];
        else  %% matrix completion
            ini = [R;zeros(n-1,1)];
        end
    end
end


