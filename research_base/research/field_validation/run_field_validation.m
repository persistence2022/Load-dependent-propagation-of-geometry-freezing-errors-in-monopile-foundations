function out=run_field_validation(mode)

if nargin==0,mode='formal';end
root=fileparts(fileparts(mfilename('fullpath')));here=fileparts(mfilename('fullpath'));addpath(fullfile(root,'pisa_model'));addpath(here);
out=fullfile(here,'runs',[mode '_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);mkdir(out);
diary(fullfile(out,'matlab_run.log'));cleanup=onCleanup(@()diary('off'));
fprintf('MATLAB %s; release %s; %s\n',version,version('-release'),computer);
G=readtable(fullfile(root,'published_data','pisa_2m_geometry.csv'),'TextType','string');
O=readtable(fullfile(root,'published_data','pisa_2m_response_metrics.csv'),'TextType','string');
writetable(G,fullfile(out,'geometry_inputs.csv'));writetable(O,fullfile(out,'published_metric_inputs.csv'));
T=[addvars(field_scenarios('Dunkirk'),repmat("Dunkirk",7,1),'Before',1,'NewVariableNames','site');...
 addvars(field_scenarios('Cowden'),repmat("Cowden",10,1),'Before',1,'NewVariableNames','site')];
writetable(T,fullfile(out,'parameter_scenarios.csv'));write_manifest(here,root,out);
check_reactions(out);curves=table;profiles=table;soilprofiles=table;metrics=table;balances=table;failures=table;
if strcmp(mode,'smoke'),piles=["DL1","CL1"];else,piles=G.pile_id';end
for id=piles
 g=G(G.pile_id==id,:);obs=O(O.pile_id==id,:);scenarios=field_scenarios(g.site);
 if strcmp(mode,'smoke'),scenarios=scenarios([1,find(scenarios.kind~="PISA",1)],:);end
 D=g.D_m;w=unique([0;logspace(-7,-2,32)'*D;linspace(.0125,.1,48)'*D;D/1000;D/10]);
 if strcmp(mode,'smoke'),w=[0;D*1e-5;D/1000;D*.01;D*.05;D/10];end
 for j=1:height(scenarios)
  sn=scenarios.scenario_id(j);[cfg,~,P,S]=field_config(id,sn,80);
  load=struct('control','head_displacement','w_m',w,'eccentricity_m',g.h_load_m);
  fprintf('START %s %s\n',id,sn);
  try
   r=solve_pisa(cfg,D,g.t_m,g.L_installed_m,load);save(fullfile(out,char(id+"_"+sn+".mat")),'r','g','P','S','load');
   v=r.response;v.pile_id=repmat(id,height(v),1);v.site=repmat(g.site,height(v),1);v.scenario_id=repmat(sn,height(v),1);curves=[curves;v];
   v=r.profile;v.pile_id=repmat(id,height(v),1);v.site=repmat(g.site,height(v),1);v.scenario_id=repmat(sn,height(v),1);profiles=[profiles;v];
   z=linspace(0,g.L_installed_m,201)';s=field_soil(z,g.site,P,S);v=struct2table(s);v.z_m=z;
   v.pile_id=repmat(id,height(v),1);v.site=repmat(g.site,height(v),1);v.scenario_id=repmat(sn,height(v),1);soilprofiles=[soilprofiles;v];
   Hinit=interp1(w,r.response.H_N,D/1000);Hult=interp1(w,r.response.H_N,D/10);K=Hinit/(D/1000)/1e6;
   metrics=[metrics;table(id,g.site,g.pair_group,sn,S.kind,K,obs.kHinit_MN_per_m,K/obs.kHinit_MN_per_m-1,...
     Hult/1e3,obs.H_0p1D_kN,Hult/(obs.H_0p1D_kN*1e3)-1,r.runtime_seconds,...
     r.domain.all_in_domain,r.domain.pisa_reaction_model,'VariableNames',{'pile_id','site','pair_group','scenario_id','kind',...
     'Kinit_pred_MN_m','Kinit_observed_MN_m','Kinit_relative_error','H01D_pred_kN','H01D_observed_kN','H01D_relative_error','runtime_s','inside_PISA_domain','PISA_reaction_model'})];
   balances=[balances;table(id,sn,max(r.response.relative_residual),r.force_balance_error_N,r.moment_balance_error_Nm,r.base_reactions(1),r.base_reactions(2),...
    'VariableNames',{'pile_id','scenario_id','max_relative_residual','force_error_N','moment_error_Nm','toe_H_N','toe_M_Nm'})];
   fprintf('DONE %s %s K %.6g MN/m H01D %.6g kN; %.3g s\n',id,sn,K,Hult/1e3,r.runtime_seconds);
  catch ME
   failures=[failures;table(id,sn,string(ME.identifier),string(ME.message),'VariableNames',{'pile_id','scenario_id','identifier','message'})];
   fprintf('FAILED %s %s: %s\n',id,sn,ME.message);
  end
  writetable(curves,fullfile(out,'predicted_curves.csv'));writetable(metrics,fullfile(out,'field_metrics.csv'));
  writetable(profiles,fullfile(out,'ultimate_profiles.csv'));writetable(soilprofiles,fullfile(out,'soil_profiles.csv'));
  writetable(balances,fullfile(out,'equilibrium_audit.csv'));if ~isempty(failures),writetable(failures,fullfile(out,'failures.csv'));end
 end
end
meta=struct('run',out,'mode',mode,'matlab',version,'release',version('-release'),'platform',computer,...
 'case_count',height(metrics),'failure_count',height(failures),'physical_piles',numel(piles),'paired_site_units',2,...
 'status','Executed unfitted comparison; published PISA-related field evidence, not blind independent PISA validation',...
 'metrics','Table3 H0.1D and secant initial stiffness atD/1000; no loading-history fit');
meta.run_relative=strrep(out(numel(root)+2:end),char(92),'/');
fid=fopen(fullfile(out,'run_metadata.json'),'w');fprintf(fid,'%s',jsonencode(meta,PrettyPrint=true));fclose(fid);
fid=fopen(fullfile(here,['current_' mode '.json']),'w');fprintf(fid,'%s',jsonencode(meta,PrettyPrint=true));fclose(fid);
fprintf('COMPLETE %s; %d successes, %d failures\n',out,height(metrics),height(failures));
end

function check_reactions(out)
rows=table;
for site=["Dunkirk","Cowden"]
 T=field_scenarios(site);S=T(find(T.kind~="PISA",1),:);z=[.2;1;4;8];D=2;soil=struct('sigma_v_eff_Pa',1e4*z,'Su_Pa',8e4+1e4*z);
 w=[-0.02;-.001;.0001;.03];[p,~,dp]=field_reaction(z,w,0*w,D,soil,S);[pm,~,~]=field_reaction(z,-w,0*w,D,soil,S);
 h=1e-9;pp=field_reaction(z,w+h,0*w,D,soil,S);pn=field_reaction(z,w-h,0*w,D,soil,S);fd=(pp-pn)/(2*h);
 rel=max(abs(dp-fd)./max(1,abs(dp)));odd=max(abs(p+pm));assert(rel<1e-5&&odd<1e-8);
 rows=[rows;table(site,rel,odd,'VariableNames',{'site','tangent_relative_error','oddness_absolute_error_N_m'})];
end
writetable(rows,fullfile(out,'reaction_module_checks.csv'));
end

function write_manifest(here,root,out)
files=[dir(fullfile(here,'*.m'));dir(fullfile(here,'PROTOCOL.txt'));dir(fullfile(root,'pisa_model','*.m'));
 dir(fullfile(root,'published_data','pisa_2m_geometry.csv'));dir(fullfile(root,'published_data','pisa_2m_response_metrics.csv'));
 dir(fullfile(root,'published_data','*soil_profiles.csv'));dir(fullfile(root,'published_data','dunkirk_interpreted_G0_profile.csv'))];
path=strings(numel(files),1);sha=path;bytes=zeros(numel(files),1);
for k=1:numel(files)
 path(k)=string(fullfile(files(k).folder,files(k).name));bytes(k)=files(k).bytes;
 fid=fopen(path(k),'r');v=fread(fid,Inf,'*uint8');fclose(fid);md=java.security.MessageDigest.getInstance('SHA-256');md.update(v);sha(k)=lower(string(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[])));
end
writetable(table(path,bytes,sha),fullfile(out,'input_code_manifest.csv'));
end
