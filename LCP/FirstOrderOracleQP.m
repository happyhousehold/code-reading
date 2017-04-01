function [fy, gy] = FirstOrderOracleQP(data,domain,y)

%%%% determine the number of nonzeros in the solution
%%%% dependent on the domain type, generate an optimal solution
tmp = data.A * y - data.b;
fy = tmp' * tmp / 2;
gy = data.A' * tmp;
