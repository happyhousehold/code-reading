% May 2008
%
% This Matlab code implements proximal gradient methods for matrix game
%  min_{x in S_n} max_{y in S_m}  y'Ax,
% where S_n is the unit simplex in R^n, and A_{ij} in [-1,1] or 0 (with probab. p).
% This code is written by Paul Tseng, Dept. Math, Univ. Washington, Seattle.
% Ref: P. Tseng, Proximal Gradient Interpolated Gradient-Projection Methods for 
% Convex-Concave Optimization: A Unified Approach, April 2008.  

for ntrials = [1];
probnum=0;

% Generate sparse matrix A
%n=input(' enter #cols n of A ');
%m=input(' enter #rows m of A ');
%p=input(' enter probability an entry of A is nonzero ');
for n = [1000 10000];
for m = [100 1000];

for p = [.01 .1];
A_flag=1;  		%0 if constant dense A;  1 if random (sparse) A
if A_flag==0
  A=([1:m]'-ones(m,1)*m/2)*[1:n] + [eye(m) zeros(m,n-m)];
  A=A/max(max(A));
else
  A=2*rand(m,n)-1;
  P=rand(m,n);
  A(P>p)=0;   	%set A_{ij}=0 if P_{ij}>p for all i,j
end
A=sparse(A);	%store A using sparse data structure

probnum=probnum+1;
%save(int2str(probnum),'A');

%eps=input(' enter termination tolerance epsilon >0 ');
for eps = [1e-3 1e-4];

if eps==1e-3 | n==1000		%skip the case of eps=1e-4 & n=10000

fprintf('ntrial= %g  n = %g  m = %g  p= %g  eps= %g \n',ntrials,n,m,p,eps);


L_flag=1;  		%0 for default constant L;  1 for adaptive L
stop_flag=1;	%0 for checking gap of y_k, v_k;  1 for checking gap of x_{k+1}, v_k  
                  %in alg 1, 2, 3.

for alg = [1 2 3 4];

%  alg = 	1 	Nesterov's 1988 method, extended by Auslender & Teboulle (1 projection)
%		2 	Nesterov's 2005 method, 1 projection variant
%		3 	Nesterov's 2005 modified method (2 projections) 
%		4 	Nemirovski's 2005 prox method
% These algorithms use h(x) = xlog(x),  grad h = log(x)+1


% Input initial feasible x.
  x=ones(n,1)/n;

  if alg==1

    k=0;
    t1=cputime;
    D1=log(n);
    D2=log(m);
    N=ceil(4*sqrt(D1*D2)/eps)-1;		%upper bound on #iterations
    L=2*D2/eps;
    if L_flag==1
      Lmax=L;
      L=L/8;
    end
    mu=eps/(2*D2);
    fprintf('Alg 1: N= %g  initial L = %g mu = %g \n',N,L,mu);

%  f_mu(x)=mu ln(sum_i e^(a_i*x/mu) /m)  where a_i is ith row of A. 

    z=x;						%initialize z to equal x
    vsum=zeros(n,1);				%initialize wted sum of dual vars to 0
    v=zeros(n,1);
    theta=1;
    gap=eps+1;
%    maxdobj=-1e30;
    while gap > eps & k <=N
      y=(1-theta)*x+theta*z;
      u=A*y;
      maxu=max(u);
      du=u-maxu;
      edu=exp(du/mu);
      sedu=sum(edu);
      dual=edu/sedu;
      grad=(dual'*A)';				%gradient of f_mu at y
      vsum=vsum+grad/theta;
      v=(1-theta)*v+theta*grad;		%new dual vars
      obj=maxu+mu*(log(sedu)-log(m));	%value of f_mu at y

      c=-grad/(theta*L)+log(max(z,1e-30));	%grad h(z) = log(z)+1
      maxc=max(c);
      dc=c-maxc;
      edc=exp(dc);

      if L_flag==1 & L<Lmax
        znew=edc/sum(edc);			%new z
        xnew=(1-theta)*x+theta*znew;	%new x
        unew=A*xnew;
        if stop_flag==1  
          u=unew;  
        end
        maxu=max(unew);
        du=unew-maxu;
        edu=exp(du/mu);
        sedu=sum(edu);
        newobj=maxu+mu*(log(sedu)-log(m));	%new obj
        Dz=log(znew)'*znew-log(z)'*z-(log(z)+1)'*(znew-z);
        if newobj <= obj+grad'*(xnew-y)+theta^2*L*Dz
          z=znew;
          x=xnew;
          theta=2*theta/(sqrt(theta^2+4)+theta);
          k=k+1;
        else
          L=min(Lmax,2*L);			%increase L & repeat iteration
        end
      else
        z=edc/sum(edc);				%new z
        x=(1-theta)*x+theta*z;		%new x

% The three theta formulas seem comparable in performance 
%      theta=theta*(sqrt(theta^2+4)-theta)/2;
        theta=2*theta/(sqrt(theta^2+4)+theta);
        k=k+1;
      end

      if (k>n) & (mod(k,5)==0)
        pobj=max(u);				%primal objective associated with y or x
%        dobj=min(grad);				%dual objective associated with y
        dobj=min(v);				%dual objective associated with y
        gap=pobj-dobj;
%        maxdobj=max(dobj,dobj);   		%This didn't help
%        gap=pobj-maxdobj;
      end  

      if (k>n) & (mod(k,1000)==0)
%        fprintf('Alg 1: k= %g L= %g obj= %g  pobj = %g dobj = %g  \n',k,L,obj,pobj,dobj);
      end

    end

    fprintf('Alg 1: k= %g L= %g obj= %g  gap = %g  pobj = %g dobj = %g  cpu= %g  k/N= %g \n',k,L,obj,gap,pobj,dobj,cputime-t1,100*k/N);

  end

  if alg==2

    k=0;
    t1=cputime;
    D1=log(n);
    D2=log(m);
    N=ceil(4*sqrt(D1*D2)/eps)-1;		%upper bound on #iterations
    L=2*D2/eps;
    if L_flag==1
      Lmax=L;
      L=L/8;
    end
    mu=eps/(2*D2);
    fprintf('Alg 2: N= %g  initial L = %g mu = %g \n',N,L,mu);

%  f_mu(x)=mu ln(sum_i e^(a_i*x/mu) /m)  where a_i is ith row of A. 

    z=ones(n,1)/n;				%initialize z to argmin_{x in S_n} h(x) 
    v=zeros(n,1);					%initialize dual vars to 0
    theta=1;
    vartheta=theta;
    gradsum=zeros(n,1);
    gap=eps+1;
    while gap > eps & k <=N
      y=(1-theta)*x+theta*z;
      u=A*y;
      maxu=max(u);
      du=u-maxu;
      edu=exp(du/mu);
      sedu=sum(edu);
      dual=edu/sedu;
      grad=(dual'*A)';				%gradient of f_mu
      v=(1-theta)*v+theta*grad;		%new dual vars
      obj=maxu+mu*(log(sedu)-log(m));	%value of f_mu
   
      gradsum=gradsum+grad/vartheta;
      c=-gradsum/L;
      maxc=max(c);
      dc=c-maxc;
      edc=exp(dc);

      if L_flag==1 & L<Lmax
        znew=edc/sum(edc);			%new z
        xnew=(1-theta)*x+theta*znew;	%new x
        unew=A*xnew;
        if stop_flag==1  
          u=unew;  
        end
        maxu=max(unew);
        du=unew-maxu;
        edu=exp(du/mu);
        sedu=sum(edu);
        newobj=maxu+mu*(log(sedu)-log(m));	%new obj
        Dz=log(znew)'*znew-log(max(z,1e-30))'*z-(log(max(z,1e-30))+1)'*(znew-z);
        if newobj <= obj+grad'*(xnew-y)+theta^2*L*Dz
          z=znew;
          x=xnew;
          theta=2*theta/(sqrt(theta^2+4)+theta);
          vartheta=theta;
          k=k+1;
%          theta=2/(k+2);
%          vartheta=2/(k+1);
        else
          L=min(Lmax,2*L);			%increase L & repeat iteration
        end
      else
        z=edc/sum(edc);				%new z
        x=(1-theta)*x+theta*z;		%new x
      
        theta=2*theta/(sqrt(theta^2+4)+theta);
         k=k+1;
%        theta=2/(k+2);
%        vartheta=2/(k+1);
      end

      if (k>n) & (mod(k,5)==0)
        pobj=max(u);				%primal objective associated with y or x
%        dobj=min(grad);				%dual objective associated with y
        dobj=min(v);				%dual objective associated with y
        gap=pobj-dobj;
      end  

      if (k>n) & (mod(k,1000)==0)
%	  fprintf('Alg 2: k= %g  L= %g  obj= %g  pobj = %g dobj = %g \n',k,L,obj,pobj,dobj);
      end

    end

    fprintf('Alg 2: k= %g  L= %g  obj= %g  gap = %g  pobj = %g dobj = %g  cpu= %g  k/N= %g\n',k,L,obj,gap,pobj,dobj,cputime-t1,100*k/N);
    
  end

  if alg==3

    k=0;
    t1=cputime;
    D1=log(n);
    D2=log(m);
    N=ceil(4*sqrt(D1*D2)/eps)-1;		%upper bound on #iterations
    L=2*D2/eps;
    if L_flag==1
      Lmax=L;
      L=L/8;
    end
    mu=eps/(2*D2);
    fprintf('Alg 3: N= %g  initial L = %g mu = %g \n',N,L,mu);

%  f_mu(x)=mu ln(sum_i e^(a_i*x/mu) /m)  where a_i is ith row of A. 

    z=ones(n,1)/n;				%initialize z to argmin_{x in S_n} h(x) 
    v=zeros(n,1);					%initialize dual vars to 0
    theta=1;
    vartheta=theta;
    gradsum=zeros(n,1);
    gap=eps+1;
    while gap > eps & k <=N
      y=(1-theta)*x+theta*z;
      u=A*y;
      maxu=max(u);
      du=u-maxu;
      edu=exp(du/mu);
      sedu=sum(edu);
      dual=edu/sedu;
      grad=(dual'*A)';				%gradient of f_mu at y
      v=(1-theta)*v+theta*grad;		%new dual vars
      obj=maxu+mu*(log(sedu)-log(m));	%value of f_mu at y 

      c=-grad/(theta*L)+log(max(z,1e-30));	%grad h(z) = log(z)+1
	maxc=max(c);
      dc=c-maxc;
      edc=exp(dc);
      hatz=edc/sum(edc);			%hat z

      gradsum=gradsum+grad/vartheta;
      c=-gradsum/L;
	maxc=max(c);
      dc=c-maxc;
      edc=exp(dc);

      if L_flag==1 & L<Lmax
        znew=edc/sum(edc);			%new z
        xnew=(1-theta)*x+theta*hatz;	%new x
        xtest=(1-theta)*x+theta*znew;
        unew=A*xnew;
        if stop_flag==1  
          u=unew;  
        end
        maxu=max(unew);
        du=unew-maxu;
        edu=exp(du/mu);
        sedu=sum(edu);
        newobj=maxu+mu*(log(sedu)-log(m));	%new x obj
        Dz=log(znew)'*znew-log(z)'*z-(log(z)+1)'*(znew-z);
        if newobj <= obj+grad'*(xtest-y)+theta^2*L*Dz
          z=znew;
          x=xnew;
          theta=2*theta/(sqrt(theta^2+4)+theta);
          vartheta=theta;
          k=k+1;
%          theta=2/(k+2);
%          vartheta=2/(k+1);
        else
          L=min(Lmax,2*L);			%increase L & repeat iteration
        end
      else
        x=(1-theta)*x+theta*hatz;		%new x
        z=edc/sum(edc);				%new z
        theta=2*theta/(sqrt(theta^2+4)+theta);
        vartheta=theta;
        k=k+1;
%        theta=2/(k+2);
%        vartheta=2/(k+1);
      end

      if (k>n) & (mod(k,5)==0)
        pobj=max(u);				%primal objective associated with y or x
%        dobj=min(grad);				%dual objective associated with y
        dobj=min(v);				%dual objective associated with y
        gap=pobj-dobj;
      end  

      if (k>n) & (mod(k,1000)==0)
%	  fprintf('Alg 3: k= %g  L= %g  obj= %g  pobj = %g dobj = %g \n',k,L,obj,pobj,dobj);
      end

    end

    fprintf('Alg 3: k= %g  L= %g  obj= %g  gap = %g  pobj = %g dobj = %g  cpu= %g  k/N= %g \n',k,L,obj,gap,pobj,dobj,cputime-t1,100*k/N);

  end

  if alg==4

    k=0;
    t1=cputime;
    D1=log(n);
    D2=log(m);
    N=ceil((D1+D2)/eps);			%upper bound on #iterations
    L=1;
    if L_flag==1
      Lmax=L;
      L=L/4;
    end
    gamma=1;
    fprintf('Alg 4: N= %g  initial L = %g gamma = %g \n',N,L,gamma);
	
%  F(x,y)= (A'*y, -A*x) with primal variables x, dual variables y. 

    xdual=ones(m,1)/m;				%initialize dual variables y
    ysum=0;
    ydualsum=0;
    gsum=0;
    gap=eps+1;
    while gap > eps

%  Extragradient iteration
      Fx=(xdual'*A)';				%gradient of phi w.r.t. x
      Fxdual=-A*x;				%gradient of phi w.r.t. y
      c=-Fx*gamma/L+log(max(x,1e-30))+1;	%grad h(x) = log(x)+1
	maxc=max(c);
      dc=c-maxc;
      edc=exp(dc);
      y=edc/sum(edc);				%y
      c=-Fxdual*gamma/L+log(max(xdual,1e-30))+1;	%grad h(x) = log(x)+1
	maxc=max(c);
      dc=c-maxc;
      edc=exp(dc);
      ydual=edc/sum(edc);			%ydual

      Fx=(ydual'*A)';				%gradient of phi w.r.t. x
      Fxdual=-A*y;				%gradient of phi w.r.t. y
      c=-Fx*gamma/L+log(max(x,1e-30))+1;	%grad h(x) = log(x)+1
	maxc=max(c);
      dc=c-maxc;
      edc=exp(dc);
      xnew=edc/sum(edc);			%new x
      c=-Fxdual*gamma/L+log(max(xdual,1e-30))+1;	%grad h(x) = log(x)+1
	maxc=max(c);
      dc=c-maxc;
      edc=exp(dc);
      xdualnew=edc/sum(edc);			%new xdual

      if L_flag==1 & L<Lmax
        Dx=log(xnew)'*xnew-log(x)'*x-(log(x)+1)'*(xnew-x);
        Dx=Dx+log(xdualnew)'*xdualnew-log(xdual)'*xdual-(log(xdual)+1)'*(xdualnew-xdual);
        if Fx'*(xnew-y)+Fxdual'*(xdualnew-ydual)+L*Dx/gamma >= 0
          x=xnew;
          xdual=xdualnew;
%  Update weighted sum of x, y
          ysum=ysum+gamma*y;
          ydualsum=ydualsum+gamma*ydual;
          gsum=gsum+gamma;
          k=k+1;
        else
          L=min(Lmax,2*L);			%increase L & repeat iteration
        end
      else
        x=xnew;
        xdual=xdualnew;
%  Update weighted sum of x, y
        ysum=ysum+gamma*y;
        ydualsum=ydualsum+gamma*ydual;
        gsum=gsum+gamma;
        k=k+1;
      end

      if (k>n) & (mod(k,5)==0)
        yavg=ysum/gsum;
        ydualavg=ydualsum/gsum;
        pobj=max(A*yavg);
        dobj=min(ydualavg'*A);
        gap=pobj-dobj;
      end  

      if (k>n) & (mod(k,500)==0)
%        fprintf('Alg 4: k= %g  L= %g  pobj = %g dobj = %g \n',k,L,pobj,dobj);
      end

    end

    fprintf('Alg 4: k= %g  L= %g  gap = %g  pobj = %g dobj = %g   cpu= %g  k/N= %g \n',k,L,gap,pobj,dobj,cputime-t1,100*k/N);

  end

%pause;

end
end
end
end 
end
end
end
