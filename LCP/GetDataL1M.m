function [data,domain,ID]=GetDataL1M;
data=[];
domain=[];
ID=[];
xdim=1024;
ydim=512;
eps=0.01;
atype=1;
xtype=1;
ntype=1;
xmax=inf;
nnul=ceil(sqrt(xdim));
nnulprob=nnul/xdim;
normalize=0;
nxny=[32,32];
while(1==1)
    strXdim=sprintf('x-dimension: %1d',xdim);
    strYdim=sprintf('y-dimension: %1d',ydim);
    strEps=sprintf('epsilon: %6.5f',eps);
    %%%
    if atype==1,
        strA='Matrix: Rand(N)';
    elseif atype==2,
        strA='Matrix: Rand(U)';
    elseif atype==3,
        strA=sprintf('Matrix: %1dx%1d grid, in-grid sensors',nxny(1),nxny(2));
    elseif atype==4,
        strA=sprintf('Matrix: %1dx%1d grid, off-grid sensors',nxny(1),nxny(2));
    elseif atype==5,
        strA=sprintf('Matrix: %1dx%1d grid, perimeter sensors',nxny(1),nxny(2));
    else
        strA=sprintf('Matrix: %1d random sources, 2*(%1d+%1d-1)=%1d perimeter sensors',xdim,nxny(1),nxny(2),...
            2*(nxny(1)+nxny(2)-1));
    end;
    if normalize==0,
        strA=sprintf('%s normalization: off',strA);
    else
        strA=sprintf('%s normalization: %1d-norm',strA,normalize);
    end;
    if xtype==1,
        strX='Signal: Rand(N)';
    elseif xtype==2,
        strX='Signal: Rand(U)';
    else
        strX='Signal: Rand(R)';
    end;
    if ntype==1,
        strN='Noise: Rand(N)';
    elseif ntype==2,
        strN='Noise: Rand(U)';
    else
        strN='Noise: Rand(R)';
    end;
    if nnul==0,
        strX=sprintf('%s nnz_prb=%6.5f inf-norm<%5.3f',strX,nnulprob,xmax);
    else
       strX=sprintf('%s nnz=%1d inf-norm<%5.3f',strX,nnul,xmax);
    end;
    %%%
    imenu=menu('Choose option',strXdim,strYdim,strEps,strA,strX,strN,'Generate','->');
    if imenu==1,
        xdim=input('x-dimension > ');
    end;
    if imenu==2,
        ydim=input('y-dimension > ');
    end;
    if imenu==3,
        eps=input('epsilon > ');
    end;
    if imenu==4,
        type=input('Type [N]ormal/[U]niform/[G]rid > ','s');
        if (type(1)=='N')|(type(1)=='n')
            atype=1;
        elseif (type(1)=='U')|(type(1)=='u')
            atype=2;
        else
            sens=input('Sensors: [i]n/[o]ut/[p]erimeter > ','s');
            if sens(1)=='i',
                atype=3;
            elseif sens(1)=='o',
                atype=4;
            else
                srcs=input('Sources randomly placed [y/n] > ','s');
                if (srcs(1)=='y')|(srcs(1)=='Y')
                    atype=6;
                else
                    atype=5;
                end;
            end;
            nxny=input('grid sizes [nx,ny] > ');
        end;
        if (atype<6)&(atype>2),
            xdim=nxny(1)*nxny(2);
        end;
        if (atype==5)|(atype==6),
            ydim=2*(sum(nxny)-1);
        end;
        normalize=input('column normalization > ');
    end;
    if imenu==5,
        type=input('Signal generation [N]ormal/[U]niform/[R]ademacher > ','s');
        if (type(1)=='N')|(type(1)=='n')
            xtype=1;
        elseif (type(1)=='U')|(type(1)=='u')
            xtype=2;
        else
            xtype=3;
        end;
        type=input('Soft sparsity control [y/n] > ','s');
        if (type(1)=='y')|(type(1)=='Y')
            nnulprob=input('nnz_prob > ');
            nnul=0;
        else
            nnul=input('nnz > ');
            nnulprob=nnul/xdim;
        end;
        xmax=input('Bound on signal entries > ');
    end;
    if imenu==6,
        type=input('Noise generation [N]ormal/[U]niform/[R]ademacher > ','s');
        if (type(1)=='N')|(type(1)=='n')
            ntype=1;
        elseif (type(1)=='U')|(type(1)=='u')
            ntype=2;
        else
            ntype=3;
        end;
    end;
    if imenu==7,
        clear pars;
        pars.eps=eps;
        pars.atype=atype;
        pars.xtype=xtype;
        pars.ntype=ntype;
        pars.normalize=normalize;
        pars.nxny=nxny;
        if nnul>0,
            pars.nnul=nnul;
        end;
        pars.nnulproba=nnulprob;
        [res,outpars]=sdatgen2(xdim,ydim,pars);
        data.xdim=size(res.A,2);
        data.ydim=size(res.A,1);
        data.b=res.b;
        data.epsilon=res.eps;
        data.sol=res.x;
        if data.xdim==data.ydim,
            data.A=res.A+0.0001*randn(1,1)*eye(data.xdim);
        else
            data.A=res.A;
        end;
        data.R=sum(abs(data.sol));
        data.genpars=outpars;
        clear res;
        data.ytype=1;
        data.multType='s';
        data.colType='s';
        %tst=input('To test matrix [y/n] > ','s');
        %if tst(1)=='y',
            %[data.s,data.Y,data.gs,data.gmx,data.trnc]=GetGammaSimpleL1M(data.A,inf,inf,'??');
        %end;
        %%%%%
        domain.xdim=data.xdim;
        domain.ydim=data.ydim;
        domain.ytype=1;
        fff=fopen('Data\ID.dat','r');
        if fff<=0,
            ID=1;
        else
            fclose(fff);
            load Data\ID.dat -ascii;
            ID=ID+1;
        end;
        save Data\ID.dat ID -ascii;
        Data.domain=domain;
        Data.data=data;
        Data.ID=ID;
        str=sprintf('save Data\\data%1d_%1d_%1d_%1d.mat Data;',atype,ydim,xdim,ID); 
        eval(str);
        
        exA = data.A' * data.A;
        Lips = max(eig(exA));
        disp(sprintf('Lips: 1.00>=%3.2f ',Lips));
        data.L = Lips;
        data.mu = 0;
        domain.n = data.xdim;
        domain.m = data.ydim;
        domain.R = data.R;
        %%%% added by Guanghui %%%
         %%% save the data file for C code
        Data.fileID=sprintf('Data%1d_%1d_%1d_%1d.ins', atype,ydim,xdim,ID);
        fInter = fopen(Data.fileID,'w');
        fprintf(fInter, 'ytype: %d\n', domain.ytype);
        fprintf(fInter, 'xdim: %d\n', data.xdim);
        fprintf(fInter, 'ydim: %d\n', data.ydim);
        fprintf(fInter, 'epsilon: %le\n', data.epsilon);
        fprintf(fInter, 'R: %le\n', data.R);
        fprintf(fInter, 'multType: %c\n', data.multType);
        fprintf(fInter, 'colType: %c\n', data.colType);
    
        fprintf(fInter, 'b:\n');
        for i=1:data.ydim,
            fprintf(fInter, '%le ', data.b(i));
        end;
        fprintf(fInter, '\nsol:\n');
        for i=1:data.xdim,
            fprintf(fInter, '%le ', data.sol(i));
        end;
        %%save A column-wise
        fprintf(fInter, '\nA:\n');
        for j=1:xdim,
            for i=1:ydim,
                fprintf(fInter, '%le ', data.A(i,j));
            end;
        end;
        fclose(fInter);
        %%% added by Guanghui %%%
  
        clear Data;
    end;        
    if imenu==8,
        break;
    end;
end;
        