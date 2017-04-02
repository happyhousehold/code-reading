function y = A_oper(x, p, m, n, k, smap)

y = zeros(m,n,k);
for ii = 1:k
    y(:,:,ii) = p.*fft2(x.*smap(:,:,ii))/sqrt(m*n);
end