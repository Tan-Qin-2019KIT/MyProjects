%%%%compute propagation Path length%%%
function [l3,l4,l5]=path1(Er,delta,xn,ym,xl,yl,ymsize)
l3=zeros(1,ymsize);l4=zeros(1,ymsize);l5=zeros(1,ymsize);
Yt=zeros(1,ymsize);b=zeros(1,ymsize);

X2=abs(xl-xn);
Yt=(-yl+ym-delta).^2;

a=-2.*X2;
b=X2.^2 +(-Yt.*Er +delta.^2) ./ (1-Er);
c=-2.*delta.^2.*X2./(1-Er);
d=(delta .*X2).^2 ./(1-Er);

lm=quartic1(a,b,c,d,ymsize);
lcoms=(1+(delta./ lm) .^2) ./Er -1;

l3=sqrt(lm.^2 +delta.^2);
l4=(ym-delta) .*sqrt(1+1./lcoms);
l5=-yl.*sqrt(1+1./lcoms);