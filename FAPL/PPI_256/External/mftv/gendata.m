clear;
% load data512;
load data256;
data.msk= p;
re= reshape(u0, 256^2,1);
data.norm= norm(re);
data.u0=u0;
data.b =zeros(256,256,8);
data.sp= sense_map;
sigma =  1e-2;
for k=1:8
    data.b(:,:,k)= data.msk .*( (fft2(sense_map(:,:,k) .* data.u0  )+ sigma * (1 + 1i) / sqrt(2) * randn(size(data.u0))));
end
save data256-2un2 data
