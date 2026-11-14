function run_design_validation(mode)

if nargin==0,mode='formal';end
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'pisa_model'));addpath(fullfile(root,'code'));
stamp=char(datetime('now','Format','yyyyMMdd_HHmmss'));out=fullfile(root,'results',['design_' mode '_' stamp]);mkdir(out);diary(fullfile(out,'matlab_run.log'));clean=onCleanup(@()diary('off'));
snapshot_execution(root,out);
fprintf('MATLAB %s; release %s; %s\n',version,version('-release'),computer);
T=readtable(fullfile(root,'sources','design_benchmarks','published_design_curves.csv'));
cases=["Clay_D1","Clay_D2","Sand_D2t_Dr55","Sand_D2t_Dr85"];
models=["Full","p_only","No_m","No_HB","No_MB","EB_limit","Exact_annulus"];
if strcmp(mode,'smoke'),cases="Sand_D2t_Dr55";models="Full";end
curves=table;metrics=table;profiles=table;balances=table;moments=table;momentmetrics=table;
for ic=1:numel(cases)
 [base,g]=design_config(cases(ic));D=g(1);t=g(2);L=g(3);h=g(4);
 w=unique([0,logspace(-7,-2,28),linspace(.0125,.1,36),.0001,.05,.1])'*D;
 if strcmp(mode,'smoke'),w=[0,.0001,.01,.1]'*D;end
 for im=1:numel(models)
  cfg=base;
  switch models(im)
   case "p_only",cfg.components=[true false false false];
   case "No_m",cfg.components(2)=false;
   case "No_HB",cfg.components(3)=false;
   case "No_MB",cfg.components(4)=false;
   case "EB_limit",cfg.shear_factor=5e3;
   case "Exact_annulus",cfg.section='exact_annulus';
  end
  loading=struct('control','head_displacement','w_m',w,'eccentricity_m',h);
  r=solve_pisa(cfg,D,t,L,loading);response=r.response;
  response.case_id=repmat(cases(ic),height(response),1);response.model=repmat(models(im),height(response),1);
  curves=[curves;response];
  rr=r.profile;rr.case_id=repmat(cases(ic),height(rr),1);rr.model=repmat(models(im),height(rr),1);profiles=[profiles;rr];
  a=pisa_assemble(cfg,D,t,L);F=zeros(a.ndof,1);F(1:2)=[r.response.H_N(end);r.response.M_theta_Nm(end)];
  [~,J,o]=pisa_residual(cfg,a,D,r.all_dofs(:,end),F);Js=a.S'*J*a.S;
  pb=sum(a.weights.*o.p_N_per_m);pm=sum(a.weights.*a.zq_m.*o.p_N_per_m);mm=sum(a.weights.*o.m_N);hb=o.Hb_N;mb=o.Mb_Nm;
  balances=[balances;table(cases(ic),models(im),pb,hb,pm,mm,L*hb,mb,r.force_balance_error_N,r.moment_balance_error_Nm,...
   max(r.response.relative_residual),norm(Js-Js','fro')/norm(Js,'fro'),r.runtime_seconds,...
   'VariableNames',{'case_id','model','distributed_H_N','base_H_N','p_resultant_moment_Nm','distributed_moment_Nm','base_H_moment_Nm','base_M_Nm','force_error_N','moment_error_Nm','residual','tangent_asymmetry','runtime_s'})];
  for regime=["ult","sd"]
   if regime=="ult",target=D/10;else,target=D/10000;end
   for reference=["3D","published_1D"]
    [x,y,hx,hy]=published_curve(T,cases(ic),reference,regime);
    limit=min(target,x(end));grid=linspace(0,limit,1001)';ref=interp1(x,y,grid,'linear');pred=interp1(w,response.H_N,grid,'pchip');
    peak=max(ref);nmae=trapz(grid,abs(pred-ref))/(limit*peak);
    ratiohalf=interp1(w,response.H_N,limit/2,'pchip')/interp1(x,y,limit/2);
    ratioend=interp1(w,response.H_N,limit,'pchip')/interp1(x,y,limit);

    delta_x=max(hx);delta_y=max(hy);low=interp1(x,y,max(grid-delta_x,0),'linear')-delta_y;
    high=interp1(x,y,min(grid+delta_x,x(end)),'linear')+delta_y;
    minimum=trapz(grid,max(max(low-pred,pred-high),0))/(limit*peak);
    maximum=trapz(grid,max(abs(pred-low),abs(pred-high)))/(limit*peak);
    row=table(cases(ic),models(im),reference,regime,limit,nmae,minimum,maximum,ratiohalf,ratioend,...
      interp1(w,response.H_N,limit,'pchip'),interp1(x,y,limit),interp1(w,response.H_N,target,'pchip'),delta_x,delta_y,...
      'VariableNames',{'case_id','model','reference','regime','domain_max_m','NMAE','NMAE_graphical_lower','NMAE_graphical_upper','ratio_half','ratio_endpoint','predicted_H_N','reference_H_N','predicted_H_at_target_N','graphical_dx_m','graphical_dH_N'});
    metrics=[metrics;row];
   end
  end
  save(fullfile(out,char(cases(ic)+"_"+models(im)+".mat")),'r','g','loading');
  fprintf('%s %s Hsd %.8g MN Hult %.8g MN time %.3g s\n',cases(ic),models(im),interp1(w,response.H_N,D/10000)*1e-6,response.H_N(end)*1e-6,r.runtime_seconds);
 end

 [xx,yy]=published_curve(T,cases(ic),'3D','ult');Href=interp1(xx,yy,D/10);
 for fraction=[.5,1]
  load=struct('H_N',linspace(0,fraction*Href,35)','M_Nm',-h*linspace(0,fraction*Href,35)');
  r=solve_pisa(base,D,t,L,load);rr=r.profile;rr.case_id=repmat(cases(ic),height(rr),1);rr.load_fraction=repmat(fraction,height(rr),1);moments=[moments;rr];
  [z,M]=published_curve(T,cases(ic),'3D','m',fraction);zg=r.profile.z_m;Mp=abs(r.profile.M_theta_Nm);grid=linspace(max(z(1),zg(1)),min(z(end),zg(end)),501)';
  ref=interp1(z,M,grid);pred=interp1(zg,Mp,grid);err=trapz(grid,abs(pred-ref))/((grid(end)-grid(1))*max(ref));
  momentmetrics=[momentmetrics;table(cases(ic),fraction,Href,err,max(Mp),max(M),...
    'VariableNames',{'case_id','load_fraction','Hult_3D_N','moment_NMAE','max_predicted_M_Nm','max_reference_M_Nm'})];
 end

 writetable(curves,fullfile(out,'predicted_curves.csv'));writetable(metrics,fullfile(out,'curve_metrics.csv'));
 writetable(profiles,fullfile(out,'ultimate_profiles.csv'));writetable(balances,fullfile(out,'reaction_balance.csv'));
 writetable(moments,fullfile(out,'matched_load_moments.csv'));writetable(momentmetrics,fullfile(out,'moment_metrics.csv'));
end
meta=struct('run',out,'mode',mode,'matlab',version,'release',version('-release'),'cases',cases,'models',models,'case_count',numel(cases),'curve_count',numel(cases)*numel(models), ...
 'claim','External published 3D numerical design cases; related PISA model family; not independent physical observations or blind test');
meta.run_relative=strrep(out(numel(root)+2:end),char(92),'/');
fid=fopen(fullfile(out,'run_metadata.json'),'w');fprintf(fid,'%s',jsonencode(meta,PrettyPrint=true));fclose(fid);
fid=fopen(fullfile(root,['current_design_' mode '.json']),'w');fprintf(fid,'%s',jsonencode(meta,PrettyPrint=true));fclose(fid);
fprintf('COMPLETE %s\n',out);
end
