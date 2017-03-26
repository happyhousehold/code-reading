%%%% This function is used to generate the instance for             %
%     min_{x \in X} \lambda_{max} (A_0 + \sum_{i=1}^{n} x_i A_i)      %
%    where X is a standard simplex in \Re^n and                     %
%    A_0, A_1, \ldots, A_n are symmetric matrices  in \Re^{m\times m}                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = NewGameInst(n, m, d)

data.A0 = sprandsym(m,d);
for i = 1:n,
    data.A{i} = sprandsym(m,d);
end;

data.exA = reshape(data.A0, [],1);
for i=1:n,
    data.exA = [data.exA,reshape(data.A{i},[],1)];
end
        
data.m =m;
data.n = n;