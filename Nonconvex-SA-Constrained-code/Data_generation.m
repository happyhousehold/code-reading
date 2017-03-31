%% Non-convex S3VM problem
%     x = ceil(sprand(data.dim, 1, data.spr));
%     y = sign(dot(x, z_sample));
%    x = ceil(sprand(data.m, data.dim, data.spr));
%x = rand(data.m, data.dim);
x1 = [ones(m,1),sprandn(m, data.dim-1,data.spr)];
y = sign(x1*data.sample);
x2 = [ones(m,1),sprandn(m, data.dim-1,data.spr)];


