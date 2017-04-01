function data = NewInst(type, m, n, R, d, s, sigma)

%%%% determine the number of nonzeros in the solution
%%%% dependent on the domain type, generate an optimal solution
data.sol = zeros(n,1);
nz = min(n, ceil(n * s));

if type == 0 %%% box intersected with simplex
    for tt = 1: nz
        pos = ceil(rand()*n);
 %       data.sol(pos) = 2*R*(rand() - 0.5);
        data.sol(pos) = rand();
    end
    
    data.sol_l1 = sum(data.sol);
    
    if data.sol_l1 > R
        data.sol = data.sol / data.sol_l1 * R;
    end
    
%    if n > 4000
%        data.A = rand(m,n);
%    else
        data.A = sprand(m,n,d);
%    end
    data.b = data.A * data.sol + sigma * randn(m,1);
end

if type == 1 %%%box
    for tt = 1: nz
        pos = ceil(rand()*n);
 %       data.sol(pos) = 2*R*(rand() - 0.5);
        data.sol(pos) = R*rand();
    end
 %   if n > 4000
 %       data.A = rand(m,n);
 %   else
        data.A = sprand(m,n,d);
 %   end
    data.b = data.A * data.sol + sigma * randn(m,1);
end

if type == 2  %%%% simplex
                
      %  x=[sort(rand(n-1,1));1];
      %  y=[0;x(1:n-1)];
      %  data.sol=R * (x-y);
      x = [zeros(n-nz,1); sort(rand(nz-1,1)); 1];
      y=[0;x(1:n-1)];
      data.sol=R * (x-y);
      
 %     if n > 4000
 %       data.A = rand(m,n);
 %     else
            data.A = sprand(m,n,d);
 %     end
      data.b = data.A * data.sol + sigma * randn(m,1);
end

if type == 3 %%% ball
     %tx = [rand(n,1)];
     for tt = 1: nz
        pos = ceil(rand()*n);
        data.sol(pos) = R*rand();
     end
     data.sol = R * data.sol / norm(data.sol);
     if n > 4000
        data.A = rand(m,n);
      else
            data.A = sprand(m,n,d);
      end
      data.b = data.A * data.sol + sigma * randn(m,1);
end

if type == 4 %%% Spectahedron for matrix completion 
    %%% first generating a random matrix with trace bounded by R
    lr = nz; %% rank
    data.MF = sprand(n, lr, d);
    data.sol = sparse(data.MF * data.MF');
    tr = trace(data.sol);
    data.sol = data.sol/tr * R; 
    
    %%% uniformly pick up m entries from X as b.
   % pos = round(rand(m,1)*n);
    
    %%% define the sparse matrix A and the right hand side b
    data.A = sparse(m, n^2);
    data.b = sparse(m,1);
%    for i = 1: m
%        posI = 1+ round(rand()*(n-1)); %% randomly pick up the entries
%        posJ = 1+ round(rand()*(posI-1));
%        data.A(i, (posI-1)*n + posJ) = 1;
%        data.A(i, posI + (posJ-1)*n) = 1;
%        data.b(i) = data.sol(posI,posJ);
%    end
  N = n^2;  
  for i = 1: m
 %   if N > 4000
 %       tmp = rand(n); 
  %  else
%        tmp = sprandsym(n,d);
%    end
 %   tmp = triu(tmp) + triu(tmp,1)';
    data.A(i,:) = reshape(sprandsym(n,d), 1, N);
  end
  data.b = data.A * reshape(data.sol,N,1) + sigma * randn(m,1);
end

