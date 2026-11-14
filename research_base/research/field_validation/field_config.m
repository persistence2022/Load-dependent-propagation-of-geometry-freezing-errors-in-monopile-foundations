function [cfg,g,P,S]=field_config(pile_id,scenario_id,n_elements)
if nargin<3,n_elements=80;end
root=fileparts(fileparts(mfilename('fullpath')));E=fullfile(root,'published_data');
G=readtable(fullfile(E,'pisa_2m_geometry.csv'),'TextType','string');g=G(G.pile_id==string(pile_id),:);assert(height(g)==1);
T=field_scenarios(g.site);S=T(T.scenario_id==string(scenario_id),:);assert(height(S)==1);
if g.site=="Dunkirk"
 T=readtable(fullfile(E,'dunkirk_interpreted_G0_profile.csv'));P=struct('Gz',T.z_m,'G',T.value);

 jj=find(diff(P.Gz)==0);assert(numel(jj)==1&&abs(P.Gz(jj)-3)<.053);P.Gz(jj:jj+1)=3;
else
 T=readtable(fullfile(E,'cowden_interpreted_soil_profiles.csv'),'TextType','string');a=T.variable=="G0_MPa";b=T.variable=="su_kPa";
 P=struct('Gz',T.z_m(a),'G',T.value(a),'Sz',T.z_m(b),'S',T.value(b));
end
cfg=pisa_config(lower(g.site));cfg.E_Pa=g.steel_E_Pa;cfg.nu=g.steel_E_Pa/(2*g.steel_G_Pa)-1;
cfg.section='exact_annulus';cfg.n_elements=n_elements;cfg.max_iterations=100;cfg.allow_extrapolation=true;
cfg.soil_function=@(z)field_soil(z,g.site,P,S);cfg.source_label='PISA field papers004/003 interpreted soil profiles; no pile-response fit';
if S.scenario_id=="PISA_p_only",cfg.components=[true false false false];end
if S.kind~="PISA"
 cfg.components=[true false false false];cfg.reaction_function=@(z,w,theta,D,soil)field_reaction(z,w,theta,D,soil,S);cfg.reaction_label=char(S.kind);
end
end
