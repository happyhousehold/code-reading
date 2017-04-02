function divU = div2(u1, u2, type)
% compute 2d finite divergence
% type = 1 for Neumann boundary condition and 2 for periodic boundary
% condition.

if nargin < 3
    type = 2;
end;

[m,n] = size(u1);

switch type
    case 1
        divU = [zeros(1,n); diff(u1,1,1)] + [zeros(m,1), diff(u2,1,2)];
    case 2
        divU = [u1(1,:)-u1(m,:); diff(u1,1,1)] + [u2(:,1)-u2(:,end), diff(u2,1,2)];
end
