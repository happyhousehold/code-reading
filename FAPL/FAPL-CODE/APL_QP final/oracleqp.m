function [f,g]=oracleqp(data,x)
 tmp=data.A * x -data.b;
 f=0.5* norm(tmp)^2;
 g= data.A' * tmp;