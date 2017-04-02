function x = At_oper(y, p, m, n, k, smap)

x = zeros(m,n);
for ii = 1:k
    x = x + (ifft2(y(:,:,ii).*p).*conj(smap(:,:,ii)));
end

x = x*sqrt(m*n);