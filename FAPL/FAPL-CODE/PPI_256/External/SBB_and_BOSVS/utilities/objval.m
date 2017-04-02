function obj = objval(U, Ux, Uy, diffs, wtv)
% calculate objective function value of U

TV = sum(sum(sqrt(Ux.*conj(Ux)+Uy.*conj(Uy))));

Fit = norm(diffs(:))^2;

obj = wtv*TV + .5*Fit;

return;    