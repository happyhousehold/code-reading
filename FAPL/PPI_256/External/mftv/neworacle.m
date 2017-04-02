 function [f,g,f0]=neworacle(data, x, control,A)
 %%% x:size(2 * data.n, 1) 
%%% g:size(2 * data.n, 1).
%%% f: f_eta (x); g: f'_eta (x); f0= f_0 (x) means eta==0.
mat_rx=reshape(x(1:data.n), data.pm, data.pn);
mat_ix=reshape(x(data.n+1: 2*data.n), data.pm, data.pn);
mat_x = mat_rx + 1i * mat_ix;

Dx = opK_TV(mat_x);
tmp = sqrt(abs(Dx(:, :, 1)).^2 + abs(Dx(:, :, 2)).^2);
Dxnorm = sum(sum(tmp));
% fq=0;
% gq=zeros(data.pm, data.pn);
% for k=1:8
%     Ax_b=data.msk .* fft2(data.sp(:,:,k) .* mat_x)-data.b(:,:,k);
%     fq= fq + norm(Ax_b,'fro')^2/data.n;
%     gq= gq+ conj(data.sp(:,:,k)) .* ifft2(data.msk .* Ax_b);
% end
Resid = A*mat_x - data.b;
fq=norm(Resid(:))^2/data.n;
gq=sum(ifft2(Resid.* data.msk) .* conj(data.sp), 3);

if data.eta == 0
f = 0.5 * fq + control.lambda * Dxnorm;    
f0=f;
ystar(:,:,1) = Dx(:,:,1)./tmp;
ystar(:,:,2) = Dx(:,:,2)./tmp;
Dtystar = opKt_TV(ystar);
reDtystar = [reshape(real(Dtystar), data.n ,1);reshape(imag(Dtystar), data.n ,1)];
gq= [reshape(real(gq),data.n,1); reshape(imag(gq), data.n,1)];
g = gq + control.lambda * reDtystar;

else
ystar = funProjy_TV(Dx / data.eta);
Dtystar = opKt_TV(ystar);
reDtystar = [reshape(real(Dtystar), data.n ,1);reshape(imag(Dtystar), data.n ,1)];
f=0.5 * fq + control.lambda * (x'*reDtystar-0.5 * data.eta * norm(reDtystar)^2);
f0 = 0.5 * fq + control.lambda * Dxnorm;


gq= [reshape(real(gq),data.n,1); reshape(imag(gq), data.n,1)];
g = gq + control.lambda * reDtystar;
end
