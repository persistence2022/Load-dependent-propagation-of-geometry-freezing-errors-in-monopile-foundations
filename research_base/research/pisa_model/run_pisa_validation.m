function run_pisa_validation()
here=fileparts(mfilename('fullpath'));out=fullfile(here,'results');if ~exist(out,'dir'),mkdir(out);end
diary(fullfile(out,'validation_matlab.log'));clean=onCleanup(@()diary('off'));rng(416,'twister');
fprintf('Actual MATLAB %s; %s; %s\n',version,computer,datetime('now'));
R=struct([]);maxcon=0;maxder=0;
for model=["cowden","dunkirk"]
 for dr=[.45,.75,.9]
  pa=pisa_parameters(model,12,10,40,dr);
  for name=["p","m","H","M"]
   a=pa.(name);x=a.xu*[0;1e-6;.01;.2;.55;.9;.999;1;1.2];[y,dy]=pisa_conic(x,a);
   q=y/a.yu;s=x/a.xu;z=x<a.xu;f=-a.n*(q-s).^2+(1-a.n)*(q-x*a.k/a.yu).*(q-1);
   maxcon=max([maxcon;abs(f(z));abs(y(1));abs(y(end)/a.yu-1);abs(pisa_conic(-x,a)+y)/a.yu]);
   h=1e-8*a.xu;j=[2:7];fd=(pisa_conic(x(j)+h,a)-pisa_conic(x(j)-h,a))/(2*h);maxder=max(maxder,max(abs(fd-dy(j)))/max(a.k,1e-12));
  end
 end
end
R=add(R,'Conic endpoints oddness and implicit residual',maxcon,1e-10);
R=add(R,'Conic analytical tangent finite difference',maxder,2e-5);

cfg=pisa_config('dunkirk');s=struct('G0_Pa',1e8,'sigma_v_eff_Pa',1e5,'Dr_fraction',.75);pa=pisa_parameters('dunkirk',10,10,40,.75);
w=.02;th=2e-6;[~,~,~,dmt,dmw]=pisa_reactions('dunkirk',w,th,10,s,pa);
h=1e-7;[~,mp]=pisa_reactions('dunkirk',w+h,th,10,s,pa);[~,mm]=pisa_reactions('dunkirk',w-h,th,10,s,pa);er1=abs((mp-mm)/(2*h)-dmw)/abs(dmw);
h=1e-10;[~,mp]=pisa_reactions('dunkirk',w,th+h,10,s,pa);[~,mm]=pisa_reactions('dunkirk',w,th-h,10,s,pa);er2=abs((mp-mm)/(2*h)-dmt)/abs(dmt);
R=add(R,'Sand dm/dw analytical cross tangent',er1,2e-5);R=add(R,'Sand dm/dtheta analytical tangent',er2,2e-5);

cfg.n_elements=1;a=pisa_assemble(cfg,10,.091,20);F=zeros(a.ndof,1);F(3)=1e6;free=3:a.ndof;q=zeros(a.ndof,1);q(free)=a.Kbeam(free,free)\F(free);
vexact=1e6*20^3/(3*a.EI)+1e6*20/a.GA;thexact=1e6*20^2/(2*a.EI);
R=add(R,'Timoshenko cantilever tip displacement',abs(q(3)/vexact-1),1e-10);R=add(R,'Timoshenko cantilever tip rotation',abs(q(4)/thexact-1),1e-10);

cfg.n_elements=12;a=pisa_assemble(cfg,10,.091,40);q=zeros(a.ndof,1);q(1:2:26)=.04*(1-a.z_m/40);q(2:2:26)=-.001;q(27:end)=.00001;
[~,J]=pisa_residual(cfg,a,10,q,zeros(a.ndof,1));v=a.S*randn(a.ndof,1);epsfd=1e-8;
rp=pisa_residual(cfg,a,10,q+epsfd*v,zeros(a.ndof,1));rm=pisa_residual(cfg,a,10,q-epsfd*v,zeros(a.ndof,1));err=norm((rp-rm)/(2*epsfd)-J*v)/norm(J*v);
R=add(R,'Assembled nonlinear directional Jacobian',err,2e-7);

G=struct([]);ns={[200,40,20,10,4,2],[120,60,24,12,6,4,3]};lens=[20,60];
href={[25.5510,25.5507,25.5491,25.5529,25.6273,25.6202],[174.3406,174.3434,174.3623,174.4155,174.7061,174.8272,175.2942]};
sref={[.5384,.5385,.5386,.5388,.5410,.5487],[.7556,.7557,.7562,.7583,.7666,.7873,.8153]};
for c=1:2
 for k=1:numel(ns{c})
  cfg=pisa_config('dunkirk');cfg.n_elements=ns{c}(k);

  if cfg.n_elements>=100,cfg.relative_tolerance=1e-6;end
  loading=struct('control','head_displacement','w_m',[0;.001;linspace(.01,1,40)'],'eccentricity_m',50);
  r=solve_pisa(cfg,10,.091,lens(c),loading);sd=r.response.H_N(2)/1e6;ult=r.response.H_N(end)/1e6;
  eH=abs(ult-href{c}(k));eS=abs(sd-sref{c}(k));balance=max(abs([r.force_balance_error_N/r.response.H_N(end),r.moment_balance_error_Nm/r.response.M_theta_Nm(end)]));
  G=[G;struct('case_id',"C"+string(1+3*(c-1)),'elements',ns{c}(k),'Hult_MN',ult,'Hult_published_MN',href{c}(k),'Hult_absolute_error_MN',eH,'Hsd_MN',sd,'Hsd_published_MN',sref{c}(k),'Hsd_absolute_error_MN',eS,'balance_relative_error',balance,'solver_relative_tolerance',cfg.relative_tolerance,'maximum_relative_residual',max(r.response.relative_residual),'runtime_seconds',r.runtime_seconds,'all_passed',eH<=.002&&eS<=.0001&&balance<=2e-7)];
  fprintf('Sand C%d, ne=%d: Hult %.8f; Hsd %.8f\n',1+3*(c-1),ns{c}(k),ult,sd);
 end
end
GT=struct2table(G);writetable(GT,fullfile(out,'published_sand_table7_verification.csv'));R=add(R,'Sand Table7 Hult maximum absolute error MN',max(GT.Hult_absolute_error_MN),.002);R=add(R,'Sand Table7 Hsd maximum absolute error MN',max(GT.Hsd_absolute_error_MN),.0001);R=add(R,'Sand global force and moment balance',max(GT.balance_relative_error),2e-7);

cfg=pisa_config('dunkirk');cfg.n_elements=20;loading=struct('H_N',[1e5;2e6],'M_Nm',[-5e6;-1e8]);ra=solve_pisa(cfg,10,.091,40,loading);cfg.soil_diameter_m=10;rb=solve_pisa(cfg,10,.091,40,loading);
R=add(R,'Frozen diameter equal actual reproduces default',norm(ra.all_dofs-rb.all_dofs,'fro')/norm(ra.all_dofs,'fro'),1e-12);
cfg.soil_diameter_m=9.5;rc=solve_pisa(cfg,10,.091,40,loading);assert(rc.domain.counterfactual_soil_diameter&&rc.beam_diameter_m==10&&rc.soil_diameter_m==9.5);

T=struct2table(R);writetable(T,fullfile(out,'implementation_verification.csv'));disp(T);
if exist(fullfile(out,'openpile_reference_curves.csv'),'file'),compare_openpile_curves();end
assert(all(T.passed),'Some numerical verification checks failed; see retained CSV/log.');
end
function R=add(R,name,value,threshold)
R=[R;struct('test',string(name),'value',value,'threshold',threshold,'passed',value<=threshold)];
end
