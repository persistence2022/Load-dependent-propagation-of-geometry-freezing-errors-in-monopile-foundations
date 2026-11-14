function report=run_delivery_checks

root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'code'),fullfile(root,'pisa_model'));
out=fullfile(root,'QA',['entrypoint_checks_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);mkdir(out);
manifest=jsondecode(fileread(fullfile(root,'pisa_model','core_manifest.json')));
names=fieldnames(manifest.files_sha256);

core={'pisa_assemble.m','pisa_beam_curvature.m','pisa_config.m','pisa_conic.m',...
 'pisa_dunkirk_soil.m','pisa_parameters.m','pisa_reactions.m','pisa_residual.m',...
 'pisa_shape.m','pisa_soil_at.m','solve_pisa.m'};
assert(numel(names)==numel(core));
for k=1:numel(core)
 key=matlab.lang.makeValidName(core{k});actual=sha256(fullfile(root,'pisa_model',core{k}));
 assert(strcmp(actual,manifest.files_sha256.(key)),'Frozen PISA core changed: %s',core{k});
end
pointers={'current_design_formal.json','current_mechanism.json',...
 fullfile('field_validation','current_formal.json'),fullfile('design_baselines','current_formal.json')};
required={{'curve_metrics.csv','predicted_curves.csv','matched_load_moments.csv'},...
 {'geometry_counterfactuals.csv','adjoint_derivatives.csv','soil_sensitivity_curves.csv'},...
 {'field_metrics.csv','numerical_checks/mesh_convergence.csv'},...
 {'curve_metrics.csv','numerical_checks/mesh_refinement.csv'}};
resolved=cell(4,1);method=cell(4,1);
for k=1:4
 [resolved{k},~,method{k}]=resolve_current_run(fullfile(root,pointers{k}),root);
 for j=1:numel(required{k}),assert(isfile(fullfile(resolved{k},required{k}{j})),'Missing required result: %s',required{k}{j});end
end
T=readtable(fullfile(resolved{2},'geometry_counterfactuals.csv'));
assert(all(ismember({'soil_parameter_domain','actual_geometry_domain','w_within_0p1D'},T.Properties.VariableNames)));
valid=false(height(T),1);
for k=1:height(T)
 [~,g]=design_config(T.case_id(k));D=T.D_m(k);
 valid(k)=D>=5&&D<=10&&g(3)/D>=2&&g(3)/D<=6&&g(4)/D>=5&&g(4)/D<=15;
end
assert(all(valid==T.actual_geometry_domain),'Stored actual-domain labels disagree with geometry.');

fake=fullfile(out,'relocation_fixture');mkdir(fake);moved=fullfile(fake,'results','fixture_run');mkdir(moved);
p=fullfile(fake,'current_test.json');writejson(p,struct('run',resolved{1},'run_relative','results/fixture_run'));
[r,~,how]=resolve_current_run(p,fake);assert(strcmp(r,char(java.io.File(moved).getCanonicalPath())));
legacy=fullfile(fake,'current_legacy.json');writejson(legacy,struct('run','Z:/unavailable/project/results/fixture_run'));
r=resolve_current_run(legacy,fake);assert(strcmp(r,char(java.io.File(moved).getCanonicalPath())));
fieldMoved=fullfile(fake,'field_validation','runs','fixture_field');mkdir(fieldMoved);
fieldPointer=fullfile(fake,'field_validation','current_formal.json');
writejson(fieldPointer,struct('run','Z:/unavailable/project/field_validation/runs/fixture_field'));
r=resolve_current_run(fieldPointer,fake);assert(strcmp(r,char(java.io.File(fieldMoved).getCanonicalPath())));

cfg=pisa_config('dunkirk');cfg.n_elements=20;
r=solve_pisa(cfg,10,.091,20,struct('control','head_displacement','w_m',[0;.001],'eccentricity_m',50));
assert(all(isfinite(r.response.H_N))&&max(r.response.relative_residual)<=cfg.relative_tolerance);
writetable(r.response,fullfile(out,'minimal_matlab_response.csv'));
report=struct('checked_at',char(datetime('now')),'matlab',version,'release',version('-release'),...
 'platform',computer,'core_files_unchanged',numel(core),'pointer_files',{pointers},...
 'resolved_runs',{resolved},'resolution_methods',{method},'relocation_checks',3,...
 'relative_path_priority',how,'geometry_rows',height(T),...
 'actual_geometry_exceedances',sum(~valid),'soil_parameter_exceedances',sum(~T.soil_parameter_domain),...
 'displacement_exceedances',sum(~T.w_within_0p1D),'minimal_solve_states',height(r.response),...
 'minimal_solve_max_relative_residual',max(r.response.relative_residual),...
 'formal_outputs_modified',false,'audit_directory',out);
writejson(fullfile(out,'checks.json'),report);disp(report);
end
function value=sha256(path)
fid=fopen(path,'r');assert(fid>=0);clean=onCleanup(@()fclose(fid));
b=fread(fid,Inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');md.update(b);
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
function writejson(path,value)
fid=fopen(path,'w');assert(fid>=0);clean=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end
