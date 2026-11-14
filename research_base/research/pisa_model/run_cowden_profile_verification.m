function run_cowden_profile_verification()
here=fileparts(mfilename('fullpath'));out=fullfile(here,'results');source=fullfile(here,'..','soil_profiles','offshore_profiles_0p05m.csv');
T=readtable(source);cfg=pisa_config('cowden');cfg.soil=table(T.z_m,T.su_kPa*1e3,T.G0_MPa*1e6,'VariableNames',{'z_m','Su_Pa','G0_Pa'});
cfg.source_label='Byrne 2020 Fig4(a,c), black offshore representative native PDF paths, soil_profiles';
ns=[200,40,20,10,4,2];hu=[12.9636,12.9639,12.9644,12.9658,12.9651,12.8759];hs=[.7395,.7396,.7396,.7399,.7420,.7492];R=struct([]);
diary(fullfile(out,'cowden_profile_matlab.log'));cl=onCleanup(@()diary('off'));
for j=1:numel(ns)
 cfg.n_elements=ns(j);cfg.relative_tolerance=1e-8;if ns(j)>=100,cfg.relative_tolerance=1e-6;end
 load=struct('control','head_displacement','w_m',[0;.001;linspace(.01,1,40)'],'eccentricity_m',50);
 r=solve_pisa(cfg,10,.091,20,load);v=[r.response.H_N(end),r.response.H_N(2)]/1e6;
 R=[R;struct('elements',ns(j),'Hult_MN',v(1),'Hult_published_MN',hu(j),'Hult_difference_pct',100*(v(1)/hu(j)-1),'Hsd_MN',v(2),'Hsd_published_MN',hs(j),'Hsd_difference_pct',100*(v(2)/hs(j)-1),'profile_reading','PDF native path interpolation; no parameter fit','solver_relative_tolerance',cfg.relative_tolerance)];
 fprintf('Cowden C1 ne=%d Hult=%.7f Hsd=%.7f MN\n',ns(j),v(1),v(2));
end
writetable(struct2table(R),fullfile(out,'cowden_table5_C1_digitized_profile.csv'));disp(struct2table(R));

cfg.n_elements=40;cfg.relative_tolerance=1e-8;E=struct([]);
for sg=[-1,0,1]
 cfg.soil.Su_Pa=max(1,(T.su_kPa+sg*T.su_graphical_bound_kPa)*1e3);

 cfg.soil.G0_Pa=max(1,(T.G0_MPa+sg*T.G0_graphical_bound_MPa)*1e6);
 r=solve_pisa(cfg,10,.091,20,load);
 E=[E;struct('graphical_bound_sign',sg,'Hult_MN',r.response.H_N(end)/1e6,'Hsd_MN',r.response.H_N(2)/1e6,'profile_sensitivity_note','su/G0 simultaneous graph-reading bounds; G0 positivity floor 1Pa; not confidence interval')];
end
writetable(struct2table(E),fullfile(out,'cowden_C1_graphical_profile_envelope.csv'));disp(struct2table(E));
fid=fopen(fullfile(out,'cowden_C4_not_run.txt'),'w');fprintf(fid,'Cowden C4 has L=60 m. The traceable Fig4 G0 profile ends at 50 m. C4 reproduction was not performed because no documented G0 extension was supplied. The original text only explicitly states linear extension of su.\n');fclose(fid);
end
