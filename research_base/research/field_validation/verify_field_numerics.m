function verify_field_numerics(formal_dir)

here=fileparts(mfilename('fullpath'));root=fileparts(here);addpath(here);addpath(fullfile(root,'pisa_model'));
addpath(fullfile(root,'code'));
if nargin==0,formal_dir=resolve_current_run(fullfile(here,'current_formal.json'),root);end
out=fullfile(formal_dir,'numerical_checks');if ~exist(out,'dir'),mkdir(out);end
diary(fullfile(out,'matlab_numerical_checks.log'));c=onCleanup(@()diary('off'));
fprintf('MATLAB %s; %s\n',version,computer);checks=table;
ids=["DL1","DL1","CL1","CL1"];sns=["PISA_full","API_phi36","PISA_full","Matlock_e005_J050"];
for j=1:numel(ids)
 for ne=[40 80 120]
  [cfg,g]=field_config(ids(j),sns(j),ne);cfg.relative_tolerance=3e-9;
  fprintf('MESH %s %s ne%d tolerance3e-9\n',ids(j),sns(j),ne);
  D=g.D_m;w=unique([0;logspace(-7,-2,32)'*D;linspace(.0125,.1,48)'*D;D/1000;D/10]);
  load=struct('control','head_displacement','w_m',w,'eccentricity_m',g.h_load_m);r=solve_pisa(cfg,D,g.t_m,g.L_installed_m,load);
  K=interp1(w,r.response.H_N,D/1000)/(D/1000)/1e6;H=r.response.H_N(end)/1e3;
  checks=[checks;table(ids(j),sns(j),ne,K,H,max(r.response.relative_residual),r.force_balance_error_N,r.moment_balance_error_Nm,...
   'VariableNames',{'pile_id','scenario_id','elements','Kinit_MN_m','H01D_kN','max_residual','force_balance_N','moment_balance_Nm'})];
 end
end
checks.Kinit_rel_to120=zeros(height(checks),1);checks.H01D_rel_to120=zeros(height(checks),1);
for j=1:numel(ids)
 ix=checks.pile_id==ids(j)&checks.scenario_id==sns(j);ref=find(ix&checks.elements==120);
 checks.Kinit_rel_to120(ix)=checks.Kinit_MN_m(ix)/checks.Kinit_MN_m(ref)-1;
 checks.H01D_rel_to120(ix)=checks.H01D_kN(ix)/checks.H01D_kN(ref)-1;
end
writetable(checks,fullfile(out,'mesh_convergence.csv'));
reg=table;
for ep=[1e-5 1e-6 1e-7]
 [cfg,g,~,S]=field_config('CL1','Matlock_e005_J050',80);cfg.relative_tolerance=3e-9;S.clay_regularization=ep;
 cfg.reaction_function=@(z,w,th,D,soil)field_reaction(z,w,th,D,soil,S);D=g.D_m;
 w=unique([0;logspace(-7,-2,32)'*D;linspace(.0125,.1,48)'*D;D/1000;D/10]);
 r=solve_pisa(cfg,D,g.t_m,g.L_installed_m,struct('control','head_displacement','w_m',w,'eccentricity_m',g.h_load_m));
 reg=[reg;table(ep,interp1(w,r.response.H_N,D/1000)/(D/1000)/1e6,r.response.H_N(end)/1e3,max(r.response.relative_residual),...
 'VariableNames',{'regularization','Kinit_MN_m','H01D_kN','max_residual'})];
end
reg.Kinit_rel_to1e7=reg.Kinit_MN_m/reg.Kinit_MN_m(end)-1;reg.H01D_rel_to1e7=reg.H01D_kN/reg.H01D_kN(end)-1;
writetable(reg,fullfile(out,'matlock_regularization.csv'));

ph=[32;36;40];ka=max((.2153*ph.^2-8.232*ph+63.657)*1000,5400);kb=max((.1978*ph.^2-10.232*ph+136.82)*1000,5400);
writetable(table(ph,ka,kb,'VariableNames',{'phi_deg','k_above_kN_m3','k_below_kN_m3'}),fullfile(out,'API_declared_k_values.csv'));
meta=struct('matlab',version,'release',version('-release'),'formal_run',formal_dir,'numerical_check_only',true,...
 'max_abs80to120_K',max(abs(checks.Kinit_rel_to120(checks.elements==80))),...
 'max_abs80to120_H',max(abs(checks.H01D_rel_to120(checks.elements==80))),...
 'Matlock_1e6_vs1e7_K',reg.Kinit_rel_to1e7(2),'Matlock_1e6_vs1e7_H',reg.H01D_rel_to1e7(2),...
 'diagnostic_tolerance',3e-9,'formal_tolerance',1e-9,...
 'tolerance_note','Initial refinement attempt stopped at1.066e-9 against1e-9; all diagnostic runs use3e-9 to avoid numerical stagnation. Formal outputs unchanged.');
fid=fopen(fullfile(out,'numerical_check_summary.json'),'w');fprintf(fid,'%s',jsonencode(meta,PrettyPrint=true));fclose(fid);disp(meta);
end
