function [x_t,er] = ProxMapping_kkt (oldx_t, Bundle, BarX, newg, newconst, prox_center, LS, control)
%%% generate the final matrix and constants then use permsolver to solve.
er=0;
if Bundle.size==0
    totalbundle = [newg; BarX.a'];
    b = [LS - newconst; BarX.b];       
else
   totalbundle = [Bundle.matrix(1 : Bundle.size, :); newg; BarX.a'];
   b = LS - Bundle.const;
   b = [b(1 : Bundle.size) ;LS - newconst; BarX.b];
end


   perm_matrix = totalbundle * totalbundle';
   perm_matrix = [perm_matrix, -eye(Bundle.size + 2)];
   perm_const = totalbundle * prox_center - b;
    
      
   [sol, solerr] = permsolver(perm_matrix, perm_const, control.M2, control);   
   if (norm(sol(1: Bundle.size+2)) >= 0 && solerr < 1e-5)
       x_t = prox_center' - sol(1: Bundle.size+2)' * totalbundle;
       x_t=x_t';
   else
%        fprintf('solerr=%d, Proxmapping_kkt is inaccurate\n', solerr);
       er=1;
       x_t = oldx_t;
   end

   