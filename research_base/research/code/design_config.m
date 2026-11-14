function [cfg,g]=design_config(caseid)

root=fileparts(fileparts(mfilename('fullpath')));
switch string(caseid)
 case "Clay_D1",g=[7.5,.068,22.5,37.5];model='cowden';
 case "Clay_D2",g=[8.75,.091,35,87.5];model='cowden';
 case "Sand_D2t_Dr55",g=[8.75,.15,35,87.5];model='dunkirk';dr=.55;
 case "Sand_D2t_Dr85",g=[8.75,.15,35,87.5];model='dunkirk';dr=.85;
 otherwise,error('Unknown design case');
end
cfg=pisa_config(model);cfg.n_elements=40;
if strcmp(model,'dunkirk')
 cfg.soil_function=@(z)pisa_dunkirk_soil(z,dr);
else
 t=readtable(fullfile(root,'soil_profiles','offshore_profiles_0p05m.csv'));
 cfg.soil=table(t.z_m,t.su_kPa*1e3,t.G0_MPa*1e6,'VariableNames',{'z_m','Su_Pa','G0_Pa'});
 cfg.source_label='Byrne 2020 Fig4 offshore black native vector profiles; graphical reading uncertainty retained';
end
end
