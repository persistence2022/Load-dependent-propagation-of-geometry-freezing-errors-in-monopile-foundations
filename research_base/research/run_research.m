function report=run_research(stage)

if nargin==0,stage='checks';end
stage=validatestring(stage,{'checks','smoke','analysis','figures','all'});
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'code'),fullfile(root,'pisa_model'),...
    fullfile(root,'field_validation'),fullfile(root,'design_baselines'));
rng(20261113,'twister');report=struct('stage',stage,'matlab',version);
switch stage
 case 'checks'
  report=run_delivery_checks;
 case 'smoke'

  run_design_validation('smoke');
  run_design_baselines('smoke');
  run_field_validation('smoke');
 case {'analysis','all'}
  run_design_validation('formal');
  run_geometry_mechanisms;
  update_domain_labels;
  run_design_baselines('formal');verify_design_baselines;
  run_field_validation('formal');verify_field_numerics;
  report=run_delivery_checks;
  if strcmp(stage,'all')
   plot_validation_figures;plot_field_figures;plot_additional_checks;
  end
 case 'figures'

  report=run_delivery_checks;
  plot_validation_figures;plot_field_figures;plot_additional_checks;
end
end
