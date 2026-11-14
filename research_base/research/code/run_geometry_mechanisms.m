function run_geometry_mechanisms

root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'pisa_model'));addpath(fullfile(root,'code'));
out=fullfile(root,'results',['mechanism_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);mkdir(out);diary(fullfile(out,'matlab_run.log'));cl=onCleanup(@()diary('off'));
snapshot_execution(root,out);
T=readtable(fullfile(root,'sources','design_benchmarks','published_design_curves.csv'));
cases=["Clay_D1","Clay_D2","Sand_D2t_Dr55","Sand_D2t_Dr85"];
geometry=table;gradients=table;depthterms=table;sensitivity=table;failures=struct('case_id',{},'task',{},'message',{});
for ic=1:numel(cases)
 [cfg,g]=design_config(cases(ic));D=g(1);t=g(2);L=g(3);h=g(4);cfg.allow_extrapolation=true;
 [x,y]=published_curve(T,cases(ic),'3D','ult');Hu=interp1(x,y,D/10);
 for fraction=[.1,.5,1]
  H=Hu*fraction;load=struct('H_N',linspace(0,H,40)','M_Nm',-h*linspace(0,H,40)');
  r=solve_pisa(cfg,D,t,L,load);a=pisa_assemble(cfg,D,t,L);q=r.all_dofs(:,end);F=zeros(a.ndof,1);F(1:2)=[H;-h*H];[~,J]=pisa_residual(cfg,a,D,q,F);
  ew=zeros(a.ndof,1);ew(1)=1;lambda=J'\ew;wrong=J\ew;
  beamD=(3*a.B'*spdiags(a.weights*a.EI,0,numel(a.weights),numel(a.weights))*a.B+a.Kshear)/D*q;
  for epsrel=[1e-3,1e-4]
   dD=D*epsrel;cp=cfg;cm=cfg;cp.soil_diameter_m=D+dD;cm.soil_diameter_m=D-dD;
   ap=pisa_assemble(cp,D,t,L);am=pisa_assemble(cm,D,t,L);
   [Rp,~,op]=pisa_residual(cp,ap,D,q,F);[Rm,~,om]=pisa_residual(cm,am,D,q,F);soilD=(Rp-Rm)/(2*dD);
   rplus=solve_pisa(cfg,D+dD,t,L,load);rminus=solve_pisa(cfg,D-dD,t,L,load);
   fd=(rplus.response.head_w_m(end)-rminus.response.head_w_m(end))/(2*dD);
   gb=-lambda'*beamD;gs=-lambda'*soilD;ga=gb+gs;gw=-wrong'*(beamD+soilD);
   gp=-sum(a.weights.*(a.N*lambda).*(op.p_N_per_m-om.p_N_per_m)/(2*dD));
   gm=-sum(a.weights.*(a.T*lambda).*(op.m_N-om.m_N)/(2*dD));
   gh=-lambda(a.base_dofs(1))*(op.Hb_N-om.Hb_N)/(2*dD);gmb=-lambda(a.base_dofs(2))*(op.Mb_Nm-om.Mb_Nm)/(2*dD);
   gradients=[gradients;table(cases(ic),fraction,H,epsrel,r.response.head_w_m(end),gb,gs,ga,fd,abs(ga-fd)/max(abs(fd),1e-12),gw,gp,gm,gh,gmb,gp+gm+gh+gmb-gs,...
    'VariableNames',{'case_id','load_fraction','H_N','relative_step','w_m','beam_derivative','soil_derivative','adjoint_derivative','central_derivative','relative_difference','incorrect_no_transpose_derivative','p_derivative','m_derivative','HB_derivative','MB_derivative','soil_partition_closure'})];
   if epsrel==1e-4
    pp=-a.weights.*(a.N*lambda).*(op.p_N_per_m-om.p_N_per_m)/(2*dD);
    mm=-a.weights.*(a.T*lambda).*(op.m_N-om.m_N)/(2*dD);
    rr=table(repmat(cases(ic),numel(pp),1),repmat(fraction,numel(pp),1),a.zq_m,pp,mm, ...
      'VariableNames',{'case_id','load_fraction','z_m','p_derivative_contribution','m_derivative_contribution'});depthterms=[depthterms;rr];
   end
  end
  for dratio=[.90,.95,1,1.05,1.10]
   for treatment=["Updated","Frozen"]
    c=cfg;if treatment=="Frozen",c.soil_diameter_m=D;end
    try
     rr=solve_pisa(c,D*dratio,t,L,load);v=rr.response.head_w_m(end);theta=rr.response.head_theta_rad(end);
     actualD=D*dratio;actualDomain=actualD>=5&&actualD<=10&&L/actualD>=2&&L/actualD<=6&&h/actualD>=5&&h/actualD<=15;
     geometry=[geometry;table(cases(ic),fraction,treatment,dratio,actualD,rr.soil_diameter_m,H,v,theta,...
       rr.domain.all_in_domain,actualDomain,v/(D*dratio)<=.1,max(rr.response.relative_residual),rr.runtime_seconds,...
       'VariableNames',{'case_id','load_fraction','treatment','D_ratio','D_m','soil_D_m','H_N','w_m','theta_rad','soil_parameter_domain','actual_geometry_domain','w_within_0p1D','residual','runtime_s'})];
    catch err
     failures(end+1)=struct('case_id',cases(ic),'task',sprintf('%s load %.2g D ratio %.2g',treatment,fraction,dratio),'message',err.message);
    end
   end
  end
 end
 for scenario=["Nominal","G0_0p8","G0_1p2","Strength_low","Strength_high"]
  c=cfg;
  if cfg.model=="cowden"
   if scenario=="G0_0p8",c.soil.G0_Pa=c.soil.G0_Pa*.8;end
   if scenario=="G0_1p2",c.soil.G0_Pa=c.soil.G0_Pa*1.2;end
   if scenario=="Strength_low",c.soil.Su_Pa=c.soil.Su_Pa*.8;end
   if scenario=="Strength_high",c.soil.Su_Pa=c.soil.Su_Pa*1.2;end
  else
   if contains(cases(ic),'55'),dr=.55;else,dr=.85;end
   gf=1;if scenario=="G0_0p8",gf=.8;end;if scenario=="G0_1p2",gf=1.2;end
   if scenario=="Strength_low",dr=dr-.05;end;if scenario=="Strength_high",dr=dr+.05;end
   c.soil_function=@(z)scaled_sand(z,dr,gf);
  end
  w=D*unique([0,logspace(-7,-2,24),linspace(.015,.1,25),.0001])';load=struct('control','head_displacement','w_m',w,'eccentricity_m',h);
  rr=solve_pisa(c,D,t,L,load);q=rr.response;q.case_id=repmat(cases(ic),height(q),1);q.scenario=repmat(scenario,height(q),1);sensitivity=[sensitivity;q];
 end
 fprintf('%s: mechanism and sensitivity analyses complete\n',cases(ic));
 writetable(geometry,fullfile(out,'geometry_counterfactuals.csv'));writetable(gradients,fullfile(out,'adjoint_derivatives.csv'));writetable(depthterms,fullfile(out,'depth_derivative_contributions.csv'));writetable(sensitivity,fullfile(out,'soil_sensitivity_curves.csv'));
end
fid=fopen(fullfile(out,'failed_cases.json'),'w');fprintf(fid,'%s',jsonencode(failures,PrettyPrint=true));fclose(fid);
meta=struct('run',out,'matlab',version,'counterfactual_rows',height(geometry),'gradient_checks',height(gradients),'failures',numel(failures), ...
 'interpretation','Deterministic model counterfactuals; no experimental causality, probabilistic confidence or PISA recalibration');
meta.run_relative=strrep(out(numel(root)+2:end),char(92),'/');
fid=fopen(fullfile(root,'current_mechanism.json'),'w');fprintf(fid,'%s',jsonencode(meta,PrettyPrint=true));fclose(fid);
fprintf('COMPLETE %s; failed cases %d\n',out,numel(failures));
end
function s=scaled_sand(z,dr,factor),s=pisa_dunkirk_soil(z,dr);s.G0_Pa=s.G0_Pa*factor;end
