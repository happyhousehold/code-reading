%Computiing gardient of cosinus smooth function

u = mvnrnd(zeros(1,data.dim),eye(data.dim))';
[obj, GR_mu] =  simulation(z+data.mu*u , z , data);
GR = GR_mu*u + 2*data.pntl*[-max(0, -z(1))+max(0, z(1)-z(2)) ; - max(0, z(1)-z(2))];