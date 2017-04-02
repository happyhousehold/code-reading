function [ux, uy] = grad2(u, type)
% comput first order forward finite difference
% type = 1 for Neumann boundary condition and 2 for periodic boundary
% condition.

if nargin < 2
    type = 2;
end;

[m,n] = size(u);

switch type
    case 1
        ux = [diff(u,1,2), zeros(m,1)];
        uy = [diff(u,1,1); zeros(1,n)];
    case 2
        ux = [diff(u,1,2), u(:,1)-u(:,end)];
        uy = [diff(u,1,1); u(1,:)-u(end,:)];
end