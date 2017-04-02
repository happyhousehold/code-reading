%%% generate the permutation matrix 
function perm_mat = perm_gen(n)
perm_mat = zeros(2^n, 2 * n);
for i = 0 : 2^n - 1
    temt=dec2bin(i)-48;
    len=length(temt);
    tem=[zeros(1,n-len),temt];
    tem_perm1=find(tem);
    tem_perm2=find(tem==0)+n;
    perm=[tem_perm1,tem_perm2];
    perm_mat(i+1,perm)=1;
end
 perm_mat = fliplr (perm_mat);
 