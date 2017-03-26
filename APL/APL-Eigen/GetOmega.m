function [omega,domega]=GetOmega(x,domain);
if (domain.type(1)=='b')|(domain.type(1)=='e'),
    domega=(x-domain.center)./domain.axes;
    omega=0.5*domega'*domega;
    domega=domega./domain.axes;
    return;
end;
if domain.type(1)=='s',
    domega=(x-domain.center)./domain.axes;
    if min(domega)<0,
        omega=inf;
        return;
    end;
    omega=sum(domega.*log(domega+1.e-100));
    domega=(log(domega+1.e-100)+1)./domain.axes;
    return;
end;
error(sprintf('unknown domain type %s',domain.type));