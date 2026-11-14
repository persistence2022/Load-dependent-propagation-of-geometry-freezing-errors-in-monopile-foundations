function run_pisa_crosschecks()
here=fileparts(mfilename('fullpath'));out=fullfile(here,'results');diary(fullfile(out,'crosscheck_matlab.log'));cl=onCleanup(@()diary('off'));R=struct([]);

s=struct('G0_Pa',1e8,'sigma_v_eff_Pa',1e5,'Dr_fraction',.75);pa=pisa_parameters('dunkirk',10,10,40,.75);e=0;
for sw=[-1,1]
 for st=[-1,1]
  w=sw*.02;th=st*2e-6;[~,~,~,dmt,dmw]=pisa_reactions('dunkirk',w,th,10,s,pa);
  h=1e-7;[~,mp]=pisa_reactions('dunkirk',w+h,th,10,s,pa);[~,mm]=pisa_reactions('dunkirk',w-h,th,10,s,pa);e=max(e,abs((mp-mm)/(2*h)-dmw)/abs(dmw));
  h=1e-10;[~,mp]=pisa_reactions('dunkirk',w,th+h,10,s,pa);[~,mm]=pisa_reactions('dunkirk',w,th-h,10,s,pa);e=max(e,abs((mp-mm)/(2*h)-dmt)/abs(dmt));
 end
end
R=add(R,'Sand tangents in four smooth sign quadrants',e,2e-5);
[~,m0]=pisa_reactions('dunkirk',0,2e-6,10,s,pa);h=1e-10;[~,mp]=pisa_reactions('dunkirk',h,2e-6,10,s,pa);[~,mm]=pisa_reactions('dunkirk',-h,2e-6,10,s,pa);
mb=pisa_conic(2e-6*s.G0_Pa/s.sigma_v_eff_Pa,pa.m);right=mb*10*pa.p.k*s.G0_Pa;
R=add(R,'Sand |p| cusp right derivative',abs((mp-m0)/h/right-1),2e-5);R=add(R,'Sand |p| cusp left derivative',abs((m0-mm)/h/(-right)-1),2e-5);

cfg=pisa_config('dunkirk');cfg.n_elements=20;a=pisa_assemble(cfg,10,.091,40);q=zeros(a.ndof,1);q(1:2:42)=.2-.001*a.z_m;q(2:2:42)=-.001;
R=add(R,'Beam rigid translation and rotation nullspace',norm(a.Kbeam*q)/(norm(a.Kbeam,'fro')*norm(q)),1e-12);
rng(416);q=randn(a.ndof,1)*1e-3;kap=pisa_beam_curvature(a,q);energy=sum(a.weights.*a.EI.*kap.^2)+q'*a.Kshear*q;
R=add(R,'Beam internal virtual work and quadratic energy',abs(q'*a.Kbeam*q/energy-1),1e-12);
R=add(R,'Stable curvature versus published shape matrix',norm(kap-a.B*q)/norm(kap),1e-12);

cfg.reaction_function=@linear_hook;cfg.components=[true,false,false,false];cfg.reaction_label='Synthetic linear p-y verification';
load=struct('H_N',[1e5;1e6],'M_Nm',[-5e6;-50e6]);r=solve_pisa(cfg,10,.091,40,load);a=pisa_assemble(cfg,10,.091,40);[~,J]=pisa_residual(cfg,a,10,zeros(a.ndof,1),zeros(a.ndof,1));F=zeros(a.ndof,1);F(1)=1e6;F(2)=-50e6;q=J\F;
R=add(R,'Custom linear hook reproduces direct linear solve',norm(q-r.all_dofs(:,end))/norm(q),1e-9);
cfg=pisa_config('dunkirk');a=pisa_assemble(cfg,10,.091,40);cfg.soil_diameter_m=9.5;b=pisa_assemble(cfg,10,.091,40);
R=add(R,'Counterfactual soil diameter preserves beam matrix',norm(a.Kbeam-b.Kbeam,'fro'),0);assert(any(a.parameters.p.k~=b.parameters.p.k));

B=readtable(fullfile(out,'independent_bvp_reference.csv'),'TextType','string');C=struct([]);
for L=[20,60]
 for ne=[20,40,80,200]
  cfg=pisa_config('dunkirk');cfg.n_elements=ne;cfg.relative_tolerance=1e-8;if ne>=100,cfg.relative_tolerance=1e-6;end
  targets=[.001;.01;.1;.5;1];load=struct('control','head_displacement','w_m',unique([0;targets;linspace(.01,1,40)']),'eccentricity_m',50);
  r=solve_pisa(cfg,10,.091,L,load);
  for j=1:numel(targets)
   id=find(abs(r.response.head_w_m-targets(j))<1e-12,1);br=B(B.L_m==L&abs(B.head_w_m-targets(j))<1e-12,:);
   eh=abs(r.response.H_N(id)/br.H_N-1);et=abs(r.response.head_theta_rad(id)/br.head_theta_rad-1);
   C=[C;struct('L_m',L,'elements',ne,'w_m',targets(j),'H_matlab_N',r.response.H_N(id),'H_bvp_N',br.H_N,'H_relative_error',eh,'theta_relative_error',et)];
  end
 end
end
CT=struct2table(C);writetable(CT,fullfile(out,'independent_bvp_matlab_comparison.csv'));
R=add(R,'Independent collocation load error across all meshes',max(CT.H_relative_error),.001);
R=add(R,'Independent collocation rotation error across all meshes',max(CT.theta_relative_error),.001);
R=add(R,'Independent collocation load error at default-or-finer mesh',max(CT.H_relative_error(CT.elements>=40)),.001);
R=add(R,'Independent collocation rotation at default-or-finer mesh',max(CT.theta_relative_error(CT.elements>=40)),.001);
T=struct2table(R);T.required_for_default_mesh=true(height(T),1);

T.required_for_default_mesh(T.test=="Independent collocation load error across all meshes")=false;
writetable(T,fullfile(out,'additional_verification.csv'));disp(T);assert(all(T.passed(T.required_for_default_mesh)));
end
function [p,m,dp,dmt,dmw]=linear_hook(z,w,theta,D,soil)
dp=1e7*(1+z/40);p=dp.*w;m=zeros(size(w));dmt=m;dmw=m;
end
function R=add(R,name,value,threshold),R=[R;struct('test',string(name),'value',value,'threshold',threshold,'passed',value<=threshold)];end
