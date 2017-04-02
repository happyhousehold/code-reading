function [Wx, Wy] = compute_w(Ux,Uy,bx,by,rho,wP)
% W = shrink2d(DU+b, tau)

UUx = Ux + bx/rho; UUy = Uy + by/rho;
V = sqrt(UUx.*conj(UUx) + UUy.*conj(UUy));
V = max(V - wP, 0) ./ max(V,eps);
Wx = V.*UUx; Wy = V.*UUy;