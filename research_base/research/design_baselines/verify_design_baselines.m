function verify_design_baselines
here=fileparts(mfilename('fullpath'));root=fileparts(here);addpath(here);addpath(fullfile(root,'pisa_model'));addpath(fullfile(root,'code'));
[~,meta]=resolve_current_run(fullfile(here,'current_formal.json'),root);out=fullfile(meta.run,'numerical_checks');mkdir(out);diary(fullfile(out,'matlab_checks.log'));c=onCleanup(@()diary('off'));
base=readtable(fullfile(meta.run,'curve_metrics.csv'),'TextType','string');checks=table;reg=table;
for id=["Clay_D1","Clay_D2","Sand_D2t_Dr55","Sand_D2t_Dr85"]
 [cfg,g,S]=baseline_config(id,2,80);cfg.relative_tolerance=1e-8;D=g(1);fprintf('REFINEMENT %s tolerance1e-8\n',id);
 w=unique([0,logspace(-7,-2,28),linspace(.0125,.1,36),.0001,.05,.1])'*D;
 r=solve_pisa(cfg,D,g(2),g(3),struct('control','head_displacement','w_m',w,'eccentricity_m',g(4)));
 for rg=["sd","ult"]
  if rg=="sd",target=D/10000;else,target=D/10;end
  val=interp1(w,r.response.H_N,target,'pchip');b=base(base.case_id==id&base.model==S.id&base.regime==rg,:);assert(height(b)==1);
  checks=[checks;table(id,S.id,rg,b.predicted_H_N,val,val/b.predicted_H_N-1,max(r.response.relative_residual),...
    'VariableNames',{'case_id','model','regime','H40_N','H80_N','relative_change','max_residual80'})];
 end
 if S.kind=="Matlock"
  [cfg,g,S]=baseline_config(id,2,40);cfg.relative_tolerance=1e-8;S.clay_regularization=1e-7;
  cfg.reaction_function=@(z,w,th,D,soil)baseline_reaction(z,w,th,D,soil,S);
  r=solve_pisa(cfg,D,g(2),g(3),struct('control','head_displacement','w_m',w,'eccentricity_m',g(4)));
  for rg=["sd","ult"]
   if rg=="sd",target=D/10000;else,target=D/10;end
   val=interp1(w,r.response.H_N,target,'pchip');b=base(base.case_id==id&base.model==S.id&base.regime==rg,:);
   reg=[reg;table(id,rg,b.predicted_H_N,val,val/b.predicted_H_N-1,'VariableNames',{'case_id','regime','H_regularization1e6_N','H_regularization1e7_N','relative_change'})];
  end
 end
end
writetable(checks,fullfile(out,'mesh_refinement.csv'));writetable(reg,fullfile(out,'matlock_regularization.csv'));
M=struct('matlab',version,'release',version('-release'),'diagnostic_tolerance',1e-8,'formal_tolerance',1e-9,...
 'max40to80',max(abs(checks.relative_change)),'max1e6to1e7',max(abs(reg.relative_change)),...
 'scope','Representative middle parameter per each4case; not all12curve mesh repeats; formal results unmodified',...
 'tolerance_note','Initial refinement attempt stagnated at3.13e-9 against3e-9. Uniform1e-8 adopted for diagnostic runs only; no physical parameter change.');
fid=fopen(fullfile(out,'verification.json'),'w');fprintf(fid,'%s',jsonencode(M,PrettyPrint=true));fclose(fid);disp(M);
end
