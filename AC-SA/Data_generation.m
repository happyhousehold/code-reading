er(index) = normrnd(0,st);
x(index,:) = rand(1,d);
y(index) = dot(x(index, :), z_sample) + er(index);