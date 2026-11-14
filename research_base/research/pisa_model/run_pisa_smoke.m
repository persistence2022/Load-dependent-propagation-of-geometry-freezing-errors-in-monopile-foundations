function run_pisa_smoke()
here=fileparts(mfilename('fullpath'));out=fullfile(here,'results');if ~exist(out,'dir'),mkdir(out);end
diary(fullfile(out,'smoke_matlab.log'));c=onCleanup(@()diary('off'));
fprintf('Actual MATLAB %s on %s\n',version,computer);
cfg=pisa_config('dunkirk');cfg.n_elements=20;
loading=struct('control','head_displacement','w_m',[0;.001;linspace(.01,1,40)'],'eccentricity_m',50);
r=solve_pisa(cfg,10,.091,20,loading);writetable(r.response,fullfile(out,'smoke_sand_C1.csv'));save(fullfile(out,'smoke_sand_C1.mat'),'r');
disp(r.response([2,end],:));fprintf('Force/moment balance residual: %.8g N, %.8g Nm\n',r.force_balance_error_N,r.moment_balance_error_Nm);
fprintf('Table7 reference Hsd=.5384 MN, Hult=25.5510 MN (200 elements).\n');
end
