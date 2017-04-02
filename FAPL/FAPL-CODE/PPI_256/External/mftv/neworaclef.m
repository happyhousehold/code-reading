 function [f,f0]=neworaclef(data, x, control,A)
%%% x:size(2 * data.n, 1) 
%%% g:size(2 * data.n, 1).
%%% f: f_eta (x); g: f'_eta (x); f0= f_0 (x) means eta==0.
mat_rx=reshape(x(1:data.n), data.pm, data.pn);
mat_ix=reshape(x(data.n+1: 2*data.n), data.pm, data.pn);
mat_x = mat_rx + 1i * mat_ix;

Dx = opK_TV(mat_x);
tmp = sqrt(abs(Dx(:, :, 1)).^2 + abs(Dx(:, :, 2)).^2);
Dxnorm = norm(tmp,1);
% fq=0;
% for k=1:8
%     Ax_b=data.msk .* fft2(data.sp(:,:,k) .* mat_x)-data.b(:,:,k);
%     fq= fq + norm(Ax_b,'fro')^2/data.n;
% end
Resid = A*mat_x - data.b;
fq=norm(Resid(:))^2/data.n;
if data.eta == 0
f = 0.5 * fq+ control.lambda * Dxnorm;    
f0=f;

else
ystar = funProjy_TV(Dx / data.eta);
f=0.5 * fq + control.lambda * (0.5/data.eta * Dxnorm^2-0.5*data.eta *norm(ystar-Dx/data.eta,'fro')^2);
f0 = 0.5*fq+ control.lambda * Dxnorm;
end