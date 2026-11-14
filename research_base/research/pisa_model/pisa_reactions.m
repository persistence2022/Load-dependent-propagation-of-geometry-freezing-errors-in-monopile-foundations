function [p,m,dp,dmt,dmw]=pisa_reactions(model,w,theta,D,soil,a)

G=soil.G0_Pa;if lower(string(model))=="cowden",stress=soil.Su_Pa;else,stress=soil.sigma_v_eff_Pa;end
vscale=stress*D./G;tscale=stress./G;
active=stress>0&G>0;vscale(~active)=1;tscale(~active)=1;
[pb,pbt]=pisa_conic(w./vscale,a.p);p=pb.*stress*D;dp=pbt.*G;
[mb,mbt]=pisa_conic(theta./tscale,a.m);
if lower(string(model))=="cowden"
 m=mb.*stress*D^2;dmt=mbt.*G*D^2;dmw=zeros(size(w));
else
 m=mb.*abs(p)*D;dmt=mbt.*abs(p)*D./tscale;dmw=mb.*sign(p).*dp*D;
end
p(~active)=0;dp(~active)=0;m(~active)=0;dmt(~active)=0;dmw(~active)=0;
end
