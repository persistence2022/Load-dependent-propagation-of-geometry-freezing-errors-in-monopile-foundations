function out=run_targeted_mesh_review(caseSubset,tolerance,baselineRoot)

if nargin<1,caseSubset=strings(0);end
if nargin<2,tolerance=1e-9;end
here=fileparts(mfilename('fullpath'));review=fileparts(here);stage=fileparts(review);
if nargin>=3
 base=char(baselineRoot);
elseif isfolder(fullfile(stage,'research','pisa_model'))
 base=fullfile(stage,'research');
else
 error('Model:MissingLocalRoot','The local research folder is required');
end
addpath(fullfile(base,'code'),fullfile(base,'pisa_model'),'-begin');
out=fullfile(review,'runs',['mesh_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);mkdir(out);
copyfile(mfilename('fullpath')+".m",fullfile(out,'entry_snapshot.m'));
diary(fullfile(out,'matlab_run.log'));cl=onCleanup(@()diary('off'));
fprintf('Targeted numerical verification, MATLAB %s, %s\n',version,computer);
[run,~]=resolve_current_run(fullfile(base,'current_mechanism.json'),base);
old=readtable(fullfile(run,'geometry_counterfactuals.csv'),'TextType','string');
curves=readtable(fullfile(base,'sources','design_benchmarks','published_design_curves.csv'));
cases=["Clay_D1","Clay_D2","Sand_D2t_Dr55","Sand_D2t_Dr85"];
if ~isempty(caseSubset),cases=cases(ismember(cases,string(caseSubset)));end
rows=table;failures=table;
for id=cases
 [cfg,g]=design_config(id);D0=g(1);D=1.1*D0;t=g(2);L=g(3);h=g(4);
 [x,y]=published_curve(curves,id,'3D','ult');H=interp1(x,y,D0/10);
 load=struct('H_N',linspace(0,H,40)','M_Nm',-h*linspace(0,H,40)');
 for treatment=["Updated","Frozen"]
  c=cfg;c.n_elements=80;c.allow_extrapolation=true;c.relative_tolerance=tolerance;
  if treatment=="Frozen",c.soil_diameter_m=D0;end
  reference=old(old.case_id==id&old.treatment==treatment&old.load_fraction==1&abs(old.D_ratio-1.1)<1e-12,:);
  assert(height(reference)==1);
  try
   r=solve_pisa(c,D,t,L,load);v=r.response.head_w_m(end);
   row=table(id,treatment,80,H,D,reference.w_m,v,v/reference.w_m-1,...
       r.response.head_theta_rad(end),max(r.response.relative_residual),...
       r.force_balance_error_N,r.moment_balance_error_Nm,r.runtime_seconds,...
       reference.actual_geometry_domain,...
       'VariableNames',{'case_id','treatment','elements','H_N','D_m','v40_m','v80_m',...
       'relative_displacement_change','theta80_rad','max_relative_residual',...
       'force_balance_N','moment_balance_Nm','runtime_s','actual_geometry_domain'});
   rows=[rows;row];
   save(fullfile(out,char(id+"_"+treatment+"_ne80.mat")),'r','g','load');
   fprintf('%s %s v40=%.12g v80=%.12g relative=%.8g\n',id,treatment,reference.w_m,v,v/reference.w_m-1);
  catch err
   failures=[failures;table(id,treatment,string(err.identifier),string(err.message),...
       'VariableNames',{'case_id','treatment','identifier','message'})];
   fprintf('FAILED %s %s: %s\n',id,treatment,err.message);
  end
  writetable(rows,fullfile(out,'response_mesh_comparison.csv'));
 end
end
bias=table;
for id=cases
 u=rows(rows.case_id==id&rows.treatment=="Updated",:);f=rows(rows.case_id==id&rows.treatment=="Frozen",:);
 if height(u)~=1||height(f)~=1,continue,end
 b40=100*(f.v40_m/u.v40_m-1);b80=100*(f.v80_m/u.v80_m-1);change=b80-b40;
 passed=abs(change)<=.5&&max(abs([u.relative_displacement_change,f.relative_displacement_change]))<=.01;
 bias=[bias;table(id,b40,b80,change,passed,u.actual_geometry_domain,...
   'VariableNames',{'case_id','bias40_percent','bias80_percent','bias_change_percentage_points','passed_mesh_screen','actual_geometry_domain'})];
end
writetable(bias,fullfile(out,'freezing_bias_mesh_comparison.csv'));
if ~isempty(failures),writetable(failures,fullfile(out,'failed_cases.csv'));end
summary=struct('matlab',version,'release',version('-release'),'platform',computer,...
 'analysis','Post-result targeted discretization audit with rules fixed before this run',...
 'cases',numel(cases),'response_states',height(rows),'complete_treatment_pairs',height(bias),'diagnostic_tolerance',tolerance,'failure_count',height(failures),...
 'passed_bias_screens',sum(bias.passed_mesh_screen),...
 'max_abs_displacement_change',max(abs(rows.relative_displacement_change)),...
 'max_abs_bias_change_percentage_points',max(abs(bias.bias_change_percentage_points)),...
 'baseline_run',run,'p10_modified',false,...
 'physical_interpretation','Same PISA model discretization check; not perturbed-geometry physical validation');
fid=fopen(fullfile(out,'summary.json'),'w');fprintf(fid,'%s',jsonencode(summary,PrettyPrint=true));fclose(fid);
fid=fopen(fullfile(review,'current_mesh_review.json'),'w');fprintf(fid,'%s',jsonencode(struct('run',out,'run_relative',strrep(out(numel(stage)+2:end),char(92),'/')),PrettyPrint=true));fclose(fid);
disp(summary);
end
