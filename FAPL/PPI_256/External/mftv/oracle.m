 function [f,g,f0]=oracle(data, x, control)
%%% x:size(2 * data.n, 1) 
%%% g:size(2 * data.n, 1).
%%% f: f_eta (x); g: f'_eta (x); f0= f_0 (x) means eta==0.
mat_rx=reshape(x(1:data.n), data.pm, data.pn);
mat_ix=reshape(x(data.n+1: 2*data.n), data.pm, data.pn);
mat_x = mat_rx + 1i * mat_ix;

Drx = opK_TV(mat_rx);
Dix = opK_TV(mat_ix);
Dx = [Drx; Dix];
tmp = sqrt(abs(Dx(:, :, 1)).^2 + abs(Dx(:, :, 2)).^2);
reDx = [reshape(Drx, 2 * data.n, 1); reshape(Dix, 2* data.n, 1)];
Dxnorm = sum(sum(tmp));
fq=0;
gq=zeros(data.pm, data.pn);
for k=1:8
    Ax_b=data.msk .* fft2(data.sp(:,:,k) .* mat_x)-data.b(:,:,k);
    reAx_b = reshape(Ax_b, data.n, 1);
    fq= fq + norm(reAx_b)^2/data.n;
    gq= gq+ conj(data.sp(:,:,k)) .* ifft2(data.msk .* Ax_b);
end

if data.eta == 0
f = 0.5 * fq + control.lambda * Dxnorm;    
f0=f;
tmp(tmp<1) = 1;
tmp(tmp==0) = 1;
ystar(:,:,1) = Dx(:,:,1)./tmp;
ystar(:,:,2) = Dx(:,:,2)./tmp;
Dtystar1 = opKt_TV(ystar(1:data.pm,:,:));
Dtystar2 = opKt_TV(ystar(data.pm+1:2*data.pm,:,:));
reDtystar = [reshape(Dtystar1, data.n ,1);reshape(Dtystar2, data.n ,1)];

gq= [reshape(real(gq),data.n,1); reshape(imag(gq), data.n,1)];
g = gq + control.lambda * reDtystar;

else
ystar = funProjy_TV(Dx / data.eta);
reystar = [reshape(ystar(1:data.pm,:,:), 2 * data.n, 1); reshape(ystar(data.pm+1 :2* data.pm,:,:), 2 * data.n, 1)];
f = 0.5 * fq + control.lambda * (reDx'* reystar - 0.5 * data.eta * norm(reystar)^2);
f0 = 0.5 * fq + control.lambda * Dxnorm;

Dtystar1 = opKt_TV(ystar(1:data.pm,:,:));
Dtystar2 = opKt_TV(ystar(data.pm+1:2*data.pm,:,:));
reDtystar = [reshape(Dtystar1, data.n ,1);reshape(Dtystar2, data.n ,1)];
gq= [reshape(real(gq),data.n,1); reshape(imag(gq), data.n,1)];
g = gq + control.lambda * reDtystar;
end



