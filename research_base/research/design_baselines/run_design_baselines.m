function out=run_design_baselines(mode)
if nargin==0,mode='formal';end
here=fileparts(mfilename('fullpath'));root=fileparts(here);addpath(here);addpath(fullfile(root,'pisa_model'));addpath(fullfile(root,'code'));addpath(fullfile(root,'field_validation'));
out=fullfile(here,'runs',[mode '_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);mkdir(out);
diary(fullfile(out,'matlab_run.log'));clean=onCleanup(@()diary('off'));
fprintf('MATLAB %s %s\n',version,computer);audit_equations(out);manifest(here,root,out);
T=readtable(fullfile(root,'sources','design_benchmarks','published_design_curves.csv'));
cases=["Clay_D1","Clay_D2","Sand_D2t_Dr55","Sand_D2t_Dr85"];curves=table;metrics=table;profiles=table;soils=table;inputs=table;balances=table;failures=table;
if strcmp(mode,'smoke'),cases=["Clay_D1","Sand_D2t_Dr55"];levels=2;else,levels=1:3;end
for id=cases
 for k=levels
  [cfg,g,S]=baseline_config(id,k);D=g(1);t=g(2);L=g(3);h=g(4);
  w=unique([0,logspace(-7,-2,28),linspace(.0125,.1,36),.0001,.05,.1])'*D;
  if strcmp(mode,'smoke'),w=[0;.0001;.001;.01;.05;.1]*D;end
  load=struct('control','head_displacement','w_m',w,'eccentricity_m',h);
  fprintf('START %s %s\n',id,S.id);
  try
   r=solve_pisa(cfg,D,t,L,load);save(fullfile(out,char(id+"_"+S.id+".mat")),'r','g','S','load');
   v=r.response;v.case_id=repmat(id,height(v),1);v.model=repmat(S.id,height(v),1);curves=[curves;v];
   v=r.profile;v.case_id=repmat(id,height(v),1);v.model=repmat(S.id,height(v),1);profiles=[profiles;v];
   z=linspace(0,L,201)';s=pisa_soil_at(cfg,z);if ~isfield(s,'Su_Pa'),s.Su_Pa=zeros(size(z));end;if ~isfield(s,'Dr_fraction'),s.Dr_fraction=zeros(size(z));end
   v=table(z,s.G0_Pa,s.Su_Pa,s.sigma_v_eff_Pa,s.Dr_fraction,'VariableNames',{'z_m','G0_Pa','Su_Pa','sigma_v_eff_Pa','Dr_fraction'});
   v.case_id=repmat(id,height(v),1);v.model=repmat(S.id,height(v),1);soils=[soils;v];
   inputs=[inputs;table(id,S.id,D,t,L,h,cfg.E_Pa,cfg.nu,string(cfg.section),cfg.shear_factor,cfg.n_elements,S.phi_deg,S.epsilon50,S.J,S.clay_regularization,...
    'VariableNames',{'case_id','model','D_m','t_m','L_m','h_m','E_Pa','nu','section','shear_factor','elements','phi_deg','epsilon50','J','regularization'})];
   balances=[balances;table(id,S.id,max(r.response.relative_residual),r.force_balance_error_N,r.moment_balance_error_Nm,r.runtime_seconds,...
    'VariableNames',{'case_id','model','max_residual','force_error_N','moment_error_Nm','runtime_s'})];
   for regime=["ult","sd"]
    if regime=="ult",target=D/10;else,target=D/10000;end
    [x,y,hx,hy]=published_curve(T,id,'3D',regime);limit=min(target,x(end));grid=linspace(0,limit,1001)';ref=interp1(x,y,grid);pred=interp1(w,r.response.H_N,grid,'pchip');
    peak=max(ref);nmae=trapz(grid,abs(pred-ref))/(limit*peak);dx=max(hx);dy=max(hy);
    low=interp1(x,y,max(grid-dx,0))-dy;high=interp1(x,y,min(grid+dx,x(end)))+dy;
    minerr=trapz(grid,max(max(low-pred,pred-high),0))/(limit*peak);maxerr=trapz(grid,max(abs(pred-low),abs(pred-high)))/(limit*peak);
    phalf=interp1(w,r.response.H_N,limit/2,'pchip');rhalf=interp1(x,y,limit/2);pend=interp1(w,r.response.H_N,limit,'pchip');rend=interp1(x,y,limit);
    metrics=[metrics;table(id,S.id,regime,limit,nmae,minerr,maxerr,phalf/rhalf,pend/rend,interp1(w,r.response.H_N,target,'pchip'),rend,dx,dy,...
      'VariableNames',{'case_id','model','regime','domain_max_m','NMAE','NMAE_graphical_lower','NMAE_graphical_upper','ratio_half','ratio_endpoint','predicted_H_N','reference_H_N','graphical_dx_m','graphical_dH_N'})];
   end
   fprintf('DONE %s %s %.6g MN H01D\n',id,S.id,r.response.H_N(end)/1e6);
  catch ME
   failures=[failures;table(id,S.id,string(ME.identifier),string(ME.message),'VariableNames',{'case_id','model','error_id','message'})];
   fprintf('FAILED %s %s %s\n',id,S.id,ME.message);
  end
  writetable(curves,fullfile(out,'predicted_curves.csv'));writetable(metrics,fullfile(out,'curve_metrics.csv'));writetable(inputs,fullfile(out,'parameter_scenarios.csv'));
  writetable(profiles,fullfile(out,'ultimate_profiles.csv'));writetable(soils,fullfile(out,'soil_profiles.csv'));writetable(balances,fullfile(out,'equilibrium_audit.csv'));
  if ~isempty(failures),writetable(failures,fullfile(out,'failures.csv'));end
 end
end
meta=struct('run',out,'mode',mode,'matlab',version,'release',version('-release'),'platform',computer,'curve_count',height(inputs),'metric_count',height(metrics),'failure_count',height(failures),...
 'claim','Unfitted API/Matlock parameter scenarios against external published3D numerical design curves; no physical-truth claim');
meta.run_relative=strrep(out(numel(root)+2:end),char(92),'/');
for p={fullfile(out,'run_metadata.json'),fullfile(here,['current_' mode '.json'])},fid=fopen(p{1},'w');fprintf(fid,'%s',jsonencode(meta,PrettyPrint=true));fclose(fid);end
fprintf('COMPLETE %s\n',out);
end
function audit_equations(out)
T=table;z=[6;8;12;16];w=[-.02;.001;-.0001;.03];soil=struct('sigma_v_eff_Pa',10090*z,'Su_Pa',1e5+1000*z);
for kind=["API","Matlock"]
 S=struct('kind',kind,'phi_deg',36,'epsilon50',.005,'J',.5,'clay_regularization',1e-6);
 [p,~,dp]=baseline_reaction(z,w,w*0,8.75,soil,S);[pf,~,df]=field_reaction(z,w,w*0,8.75,soil,S);
 eq=max(abs(p-pf));de=max(abs(dp-df));hh=1e-9;fd=(baseline_reaction(z,w+hh,0*w,8.75,soil,S)-baseline_reaction(z,w-hh,0*w,8.75,soil,S))/(2*hh);
 rel=max(abs(fd-dp)./max(1,abs(dp)));assert(eq<1e-7&&de<1e-5&&rel<1e-5);
 T=[T;table(kind,eq,de,rel,'VariableNames',{'kind','formula_difference_N_m','tangent_difference','finite_difference_relative_error'})];
end
writetable(T,fullfile(out,'reaction_module_checks.csv'));
end
function manifest(here,root,out)
f=[dir(fullfile(here,'*.m'));dir(fullfile(here,'PROTOCOL.txt'));dir(fullfile(root,'pisa_model','*.m'));dir(fullfile(root,'code','design_config.m'));dir(fullfile(root,'code','published_curve.m'));
 dir(fullfile(root,'field_validation','field_reaction.m'));dir(fullfile(root,'sources','design_benchmarks','published_design_curves.csv'));dir(fullfile(root,'soil_profiles','offshore_profiles_0p05m.csv'))];
p=strings(numel(f),1);sha=p;b=zeros(numel(f),1);
for j=1:numel(f),p(j)=string(fullfile(f(j).folder,f(j).name));b(j)=f(j).bytes;fid=fopen(p(j),'r');v=fread(fid,Inf,'*uint8');fclose(fid);md=java.security.MessageDigest.getInstance('SHA-256');md.update(v);sha(j)=lower(string(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[])));end
writetable(table(p,b,sha,'VariableNames',{'file','bytes','sha256'}),fullfile(out,'input_code_manifest.csv'));
end
