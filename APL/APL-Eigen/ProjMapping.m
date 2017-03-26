function [x,omeganew,domeganew]=ProjMapping(xi,domain);
% computes x=argmin [\langle \xi,x\rangle + \omega(x)]
% and omeganew=\omega(x), domeganew=\omega'(x)
domeganew=xi.*domain.axes;
if (domain.type(1)=='e')|(domain.type(1)=='b')
    if domain.type(1)=='e',
        if domain.type(2)=='p',
            domeganew=max(0,domeganew);
        end;
        nrm=max(norm(domeganew),1);
        x=domain.center-(domeganew/nrm)*domain.axes;
    end;
    if domain.type(1)=='b',
        if domain.type(2)=='p',
            x=(domeganew<0).*min(1,-domeganew);
        else
            x=(domeganew<0).*min(1,-domeganew)-(domeganew>0).*min(1,domeganew);
        end;
        x=domain.center+x.*domain.axes;
    end;
end;
if domain.type(1)=='s',
    mx=max(-domeganew-1);
    x=exp(-domeganew-1-mx);
    if domain.type(2)=='e',
        x=x/sum(x);
    elseif log(sum(x))+mx>0,
        x=x/sum(x);
    else
        x=x*exp(mx);
    end;
    x=domain.center+x.*domain.axes;
end;
[omeganew,domeganew]=GetOmega(x,domain);
