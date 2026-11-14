function [p,dp,meta]=soil_py(cfg,z,y,D)

arguments
 cfg struct
 z (:,1) double {mustBeNonnegative,mustBeFinite}
 y (:,1) double {mustBeFinite}
 D (1,1) double {mustBePositive,mustBeFinite}
end
assert(numel(z)==numel(y));
if strcmp(cfg.soil_model,'linear')
 p=cfg.linear_k_N_m2*y;dp=ones(size(y))*cfg.linear_k_N_m2;meta=struct('capacity',inf(size(y)),'layer',zeros(size(y)));return
end
tops=[0 cfg.layer_bottom_m(1:end-1)];n=numel(z);p=zeros(n,1);dp=p;cap=p;idx=p;
for k=1:numel(tops)
 mask=z>=tops(k)&z<=cfg.layer_bottom_m(k);
 if k>1,mask=mask&z>tops(k);end
 if ~any(mask),continue;end
 zz=z(mask);yy=y(mask);sig=zeros(size(zz));
 for j=1:numel(tops)
  sig=sig+cfg.gamma_N_m3(j).*max(0,min(zz,cfg.layer_bottom_m(j))-tops(j));
 end
 if cfg.layer_type(k)==1
  cu=cfg.cu_Pa(k);pu=min((3*cu+sig)*D+0.5*cu*zz,9*cu*D);
  yc=2.5*cfg.e50(k)*D;r=yy/yc;ep=cfg.clay_regularization;
  base=0.5*pu.*r./(r.^2+ep^2).^(1/3);
  tang=0.5*pu/yc.*(ep^2+r.^2/3)./(r.^2+ep^2).^(4/3);
  sat=abs(base)>=pu;base(sat)=sign(base(sat)).*pu(sat);tang(sat)=0;
  if cfg.py_linear
   ff=cfg.linearization_fraction;klin=pu./(8*yc*ff^2);base=klin.*yy;tang=klin;
  end
  p(mask)=base;dp(mask)=tang;cap(mask)=pu;
 else
  ph=cfg.phi_deg(k)*pi/180;al=ph/2;be=pi/4+ph/2;ko=0.4;ka=tan(pi/4-ph/2)^2;
  pu1=sig.*(ko*zz*tan(ph)*sin(be)/(tan(be-ph)*cos(al)) + tan(be)/tan(be-ph)*(D+zz*tan(be)*tan(al)) + ko*zz*tan(be)*(tan(ph)*sin(be)-tan(al))-ka*D);
  pu2=D*sig.*(ka*(tan(be)^8-1)+ko*tan(ph)*tan(be)^4);
  A=max(0.9,3-0.8*zz/D);pu=A.*min(pu1,pu2);kz=cfg.k_N_m3(k)*zz;
  aa=kz.*yy./max(pu,realmin);th=tanh(aa);p(mask)=pu.*th;dp(mask)=kz.*(1-th.^2);cap(mask)=pu;
  if cfg.py_linear
   ff=cfg.linearization_fraction;klin=(ff/atanh(ff))*kz;p(mask)=klin.*yy;dp(mask)=klin;
  end
 end
 idx(mask)=k;
end
assert(all(idx>0),'physics:SoilCoverage','Depth outside profile.');
meta=struct('capacity',cap,'layer',idx);
end
