function [cfg,g,S]=baseline_config(case_id,k,n_elements)
if nargin<3,n_elements=40;end
[cfg,g]=design_config(case_id);cfg.n_elements=n_elements;cfg.max_iterations=100;
S=struct('kind',"",'phi_deg',nan,'epsilon50',nan,'J',.5,'clay_regularization',1e-6,'id',"");
if cfg.model=="dunkirk"
 ph=[32 36 40];S.kind="API";S.phi_deg=ph(k);S.id="API_phi"+ph(k);
else
 ep=[.004 .005 .007];S.kind="Matlock";S.epsilon50=ep(k);S.id="Matlock_e"+compose('%03d',round(ep(k)*1000))+"_J050";
 cfg.soil.sigma_v_eff_Pa=(21.19-9.81)*1e3*cfg.soil.z_m;
end
cfg.components=[true false false false];cfg.reaction_function=@(z,w,t,D,soil)baseline_reaction(z,w,t,D,soil,S);
cfg.reaction_label=char(S.id);
end
