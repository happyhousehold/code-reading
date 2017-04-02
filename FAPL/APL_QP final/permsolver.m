function [sol, solerr] = permsolver(perm_matrix, perm_const, limit_matrix, control)
%%% this function is to solve problem: perm_matrix * x = perm_const
%%% such that x_i > = 0 and half of x_i ==0, size(perm_matrix)=(n, 2*n),
%%% x: size(2*n, 1), perm_const; size(2*n, 1). sol: size(2*n, 1).
warning off;
n = size(perm_matrix, 1);
if n == control.bundle_limit + 2;
    perm_mat= limit_matrix;
else
perm_mat = perm_gen(n);
end

sol = zeros(2 * n, 1);
solerr=1;
tmpsol=zeros(2 * n, 1);
tmperr=1;
for k=1:2^n;
    ch=find(perm_mat(k,:));
    nch=find(perm_mat(k,:)==0); 
    perm_A = perm_matrix(: , ch);
  if norm(perm_A)>0
        perm_sol=perm_A\perm_const;
  end
  
  if min(perm_sol)>=0
      tmpsol(ch)= perm_sol;
      tmpsol(nch)=0;
      tmperr = norm(abs(perm_matrix * tmpsol - perm_const));
      if tmperr< 1e-12
          sol = tmpsol;
          solerr = tmperr;
          break;
      else
          if tmperr < solerr
              sol=tmpsol;
              solerr = tmperr;
          end
      end
  end
end
warning on;

      



    