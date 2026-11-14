function [Nv,Nt,Bk]=pisa_shape(s,h)

N=[1-3*s.^2+2*s.^3,h*(s-2*s.^2+s.^3),3*s.^2-2*s.^3,h*(-s.^2+s.^3)];
dN=[(-6*s+6*s.^2)/h,1-4*s+3*s.^2,(6*s-6*s.^2)/h,-2*s+3*s.^2];
ddN=[(-6+12*s)/h^2,(-4+6*s)/h,(6-12*s)/h^2,(-2+6*s)/h];
Nv=[N(1),N(2),N(2)+N(4),N(3),N(4)];
Nt=[dN(1),dN(2),dN(2)+dN(4)-1,dN(3),dN(4)];
Bk=[ddN(1),ddN(2),ddN(2)+ddN(4),ddN(3),ddN(4)];
end
