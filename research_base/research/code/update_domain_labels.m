function update_domain_labels

root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'code'));[~,s]=resolve_current_run(fullfile(root,'current_mechanism.json'),root);f=fullfile(s.run,'geometry_counterfactuals.csv');T=readtable(f);
if ismember('geometry_in_calibration_domain',T.Properties.VariableNames)
 T=renamevars(T,'geometry_in_calibration_domain','soil_parameter_domain');
end
valid=false(height(T),1);
for k=1:height(T),[~,g]=design_config(T.case_id(k));D=T.D_m(k);valid(k)=D>=5&&D<=10&&g(3)/D>=2&&g(3)/D<=6&&g(4)/D>=5&&g(4)/D<=15;end
T.actual_geometry_domain=valid;writetable(T,f);
record=struct('soil_parameter_domain_exceedances',sum(~T.soil_parameter_domain),'actual_geometry_exceedances',sum(~valid),'displacement_range_exceedances',sum(~T.w_within_0p1D));
fid=fopen(fullfile(s.run,'domain_summary.json'),'w');fprintf(fid,'%s',jsonencode(record,PrettyPrint=true));fclose(fid);disp(record);
end
