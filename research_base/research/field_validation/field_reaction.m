function [p,m,dp,dmt,dmw]=field_reaction(z,w,theta,D,soil,S)

z=z(:);w=w(:);sig=soil.sigma_v_eff_Pa;
if S.kind=="API"
 ph=S.phi_deg*pi/180;al=ph/2;be=pi/4+ph/2;ko=.4;ka=tan(pi/4-ph/2)^2;
 c1=ko*tan(ph)*sin(be)/(tan(be-ph)*cos(al))+tan(be)/tan(be-ph)*tan(be)*tan(al)+ko*tan(be)*(tan(ph)*sin(be)-tan(al));
 c2=tan(be)/tan(be-ph)-ka;c3=ka*(tan(be)^8-1)+ko*tan(ph)*tan(be)^4;
 A=max(.9,3-.8*z/D);cap=A.*min(sig.*(c1*z+c2*D),c3*sig*D);
 ph=S.phi_deg;kbelow=max((.1978*ph^2-10.232*ph+136.82)*1000,5400)*1000;
 kabove=max((.2153*ph^2-8.232*ph+63.657)*1000,5400)*1000;
 kz=(kabove*(z<5.4)+kbelow*(z>=5.4)).*z;
 a=kz.*w./max(cap,realmin);th=tanh(a);p=cap.*th;dp=kz.*(1-th.^2);
else
 assert(S.kind=="Matlock");su=soil.Su_Pa;cap=min((3*su+sig)*D+S.J*su.*z,9*su*D);
 y50=2.5*S.epsilon50*D;r=w/y50;ep=S.clay_regularization;
 p=.5*cap.*r./(r.^2+ep^2).^(1/3);dp=.5*cap/y50.*(ep^2+r.^2/3)./(r.^2+ep^2).^(4/3);
 sat=abs(p)>=cap;p(sat)=sign(p(sat)).*cap(sat);dp(sat)=0;
end
m=zeros(size(z));dmt=m;dmw=m;
end
