stream1 = RandStream('mt19937ar','Seed',seed_val);
RandStream.setDefaultStream(stream1);
r1 = randi(stream1,N_iter,1,num_rnd);
%r1 = rand(stream1, 1, num_rnd);