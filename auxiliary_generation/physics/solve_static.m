function r=solve_static(cfg,D,t,L,mode)

arguments
 cfg struct
 D (1,1) double {mustBePositive,mustBeFinite}
 t (1,1) double {mustBePositive,mustBeFinite}
 L (1,1) double {mustBePositive,mustBeFinite}
 mode {mustBeMember(mode,{'updated','frozen','updated_linear','frozen_linear'})}
end
tic;assert(2*t<D&&L<=cfg.layer_bottom_m(end));assert(cfg.max_element_m>0);
I=pi/64*(D^4-(D-2*t)^4);EI=cfg.E_Pa*I;
Ds=D;if startsWith(mode,'frozen'),Ds=cfg.D0_m;end
cfg.py_linear=endsWith(mode,'_linear');
if cfg.py_linear,cfg.load_steps=1;end

z=linspace(0,L,ceil(L/cfg.max_element_m)+1)';ne=numel(z)-1;nd=2*numel(z);
gp=[-sqrt(3/5),0,sqrt(3/5)];gw=[5/9,8/9,5/9];
zg=[];weights=[];ge=[];
for e=1:ne
 cuts=unique([z(e),cfg.layer_bottom_m(cfg.layer_bottom_m>z(e)&cfg.layer_bottom_m<z(e+1)),z(e+1)]);
 for q=1:numel(cuts)-1
  dh=cuts(q+1)-cuts(q);zg=[zg;cuts(q)+(gp'+1)*dh/2];weights=[weights;gw'*dh/2];ge=[ge;repmat(e,3,1)];
 end
end
ng=numel(zg);K=sparse(nd,nd);B=sparse(ng,nd);Bs=sparse(2*ne,nd);zs=zeros(2*ne,1);
for e=1:ne
 h=z(e+1)-z(e);d=2*e-1:2*e+2;
 ke=EI/h^3*[12 6*h -12 6*h;6*h 4*h^2 -6*h 2*h^2;-12 -6*h 12 -6*h;6*h 2*h^2 -6*h 4*h^2];K(d,d)=K(d,d)+ke;
 for j=find(ge==e)'
  s=(zg(j)-z(e))/h;B(j,d)=[1-3*s^2+2*s^3,h*(s-2*s^2+s^3),3*s^2-2*s^3,h*(-s^2+s^3)];
 end
 for a=1:2
  s=a-1;j=2*(e-1)+a;Bs(j,d)=[(-6+12*s)/h^2,(-4+6*s)/h,(6-12*s)/h^2,(-2+6*s)/h];zs(j)=z(e)+s*h;
 end
end
f=zeros(nd,1);f(1:2)=[cfg.H_N;cfg.M_Nm];u=zeros(nd,1);allok=true;iters=0;res=inf;
if strcmp(cfg.tip_boundary,'fixed'),free=1:nd-2;else,free=1:nd;end
scf=ones(nd,1)*max(abs(cfg.H_N),abs(cfg.M_Nm)/max(L,1));scf(2:2:end)=scf(2:2:end)*L;scf=max(scf,1);
for step=1:cfg.load_steps
 target=f*(step/cfg.load_steps);success=false;
 for it=1:cfg.max_iterations
  [p,dp]=soil_py(cfg,zg,B*u,Ds);fs=B'*(weights.*p);rr=target-K*u-fs;res=max(abs(rr(free))./scf(free));iters=iters+1;
  if res<cfg.newton_tol,success=true;break;end
  KT=K+B'*spdiags(weights.*dp,0,ng,ng)*B;KK=KT(free,free);ss=1./sqrt(diag(KK));S=spdiags(ss,0,numel(free),numel(free));delta=zeros(nd,1);delta(free)=ss.*((S*KK*S)\(ss.*rr(free)));
  baseline=norm(rr(free)./scf(free));a=1;
  for ls=1:20
   trial=u+a*delta;pt=soil_py(cfg,zg,B*trial,Ds);rt=target-K*trial-B'*(weights.*pt);
   if norm(rt(free)./scf(free))<baseline,break;end
   a=a/2;
  end
  u=trial;
 end
 if ~success,allok=false;break;end
end
[p,dp,meta]=soil_py(cfg,zg,B*u,Ds);fs=B'*(weights.*p);KT=K+B'*spdiags(weights.*dp,0,ng,ng)*B;
head=1:2;inner=setdiff(free,head);kh=full(KT(head,head)-KT(head,inner)*(KT(inner,inner)\KT(inner,head)));
moment=EI*(Bs*u);stress=abs(moment)*D/(2*I);force_soil=sum(weights.*p);moment_soil=sum(weights.*p.*zg);

balance=max(abs([force_soil-cfg.H_N,moment_soil-cfg.M_Nm])./[max(abs(cfg.H_N),1),max(abs(cfg.M_Nm),1)]);
r=struct('D_m',D,'t_m',t,'L_m',L,'mode',mode,'soil_D_m',Ds,'rotation_rad',u(2),...
 'displacement_m',u(1),'max_stress_Pa',max(stress),'mass_kg',cfg.rho_kg_m3*pi/4*(D^2-(D-2*t)^2)*L,...
 'converged',allok,'residual_relative',res,'equilibrium_relative',balance,'iterations',iters,...
 'elements',ne,'min_element_m',min(diff(z)),'max_element_m',max(diff(z)),...
 'integration_points',ng,'head_tangent_stiffness_2x2',kh,'elapsed_seconds',toc);
if cfg.return_profile
 r.depth_m=z;r.nodal_displacement_m=u(1:2:end);r.nodal_rotation_rad=u(2:2:end);r.stress_depth_m=zs;r.bending_moment_Nm=moment;r.gauss_depth_m=zg;r.soil_reaction_N_m=p;r.soil_capacity_N_m=meta.capacity;r.soil_displacement_m=B*u;
end
if ~allok,error('physics:Nonconvergence','Newton did not converge: D=%g t=%g L=%g residual=%g',D,t,L,res);end
end
