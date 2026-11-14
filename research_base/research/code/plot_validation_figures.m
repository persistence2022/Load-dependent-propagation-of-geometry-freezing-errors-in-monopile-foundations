function plot_validation_figures(selected)

if nargin==0,selected=strings(0);else,selected=string(selected);end
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'code'),fullfile(root,'pisa_model'));
out=fullfile(root,'figures');if ~isfolder(out),mkdir(out);end
[~,d]=resolve_current_run(fullfile(root,'current_design_formal.json'),root);[~,m]=resolve_current_run(fullfile(root,'current_mechanism.json'),root);
C=readtable(fullfile(d.run,'predicted_curves.csv'));E=readtable(fullfile(d.run,'curve_metrics.csv'));P=readtable(fullfile(d.run,'matched_load_moments.csv'));
B=readtable(fullfile(d.run,'reaction_balance.csv'));G=readtable(fullfile(m.run,'geometry_counterfactuals.csv'));A=readtable(fullfile(m.run,'adjoint_derivatives.csv'));
Z=readtable(fullfile(m.run,'depth_derivative_contributions.csv'));S=readtable(fullfile(m.run,'soil_sensitivity_curves.csv'));
T=readtable(fullfile(root,'sources','design_benchmarks','published_design_curves.csv'));
cases=["Clay_D1","Clay_D2","Sand_D2t_Dr55","Sand_D2t_Dr85"];
labels={'Clay D1','Clay D2','Sand 55%','Sand 85%'};
blue=[0 114 178]/255;orange=[213 94 0]/255;green=[0 158 115]/255;purple=[204 121 167]/255;gray=[.48 .5 .52];
cols=[blue;orange;green;purple;[.34 .7 .9];[.9 .65 0];[.15 .15 .15]];
styles={'-','--','-.',':'};reg={};
set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman',...
 'defaultAxesFontSize',9,'defaultTextFontSize',9,'defaultAxesLineWidth',.7,'defaultLineLineWidth',1.25,...
 'defaultAxesTickDir','out','defaultAxesBox','off','defaultLegendBox','off','defaultLegendFontSize',8);

f=canvas(2,2,18,13);
geom=readtable(fullfile(root,'published_data','pisa_2m_geometry.csv'));
nexttile;hold on;scatter([7.5 8.75],[3 4],55,blue,'o','filled');scatter([8.75 8.75],[4 4],55,orange,'s');
scatter([2 2],[10.35 10.6]/2,45,green,'^','filled');scatter([2 2],[10.61 10.57]/2,45,purple,'d','filled');
xlabel('Diameter, {\it D} (m)');ylabel('Embedment ratio, {\it L/D}');xlim([1 10]);ylim([1.8 5.8]);
legend({'Clay design','Sand design','Cowden pair','Dunkirk pair'},'Location','southwest');panel(1);
nexttile;hold on;scatter([7.5 8.75],[.068 .091],55,blue,'o','filled');scatter(8.75,.15,55,orange,'s','filled');scatter(2,.025,55,green,'^','filled');scatter(2,.038,55,purple,'d','filled');
xlabel('Diameter, {\it D} (m)');ylabel('Wall thickness, {\it t} (m)');xlim([1 10]);ylim([0 .17]);panel(2);
nexttile;bar([4 4 2],'FaceColor',blue);xticks(1:3);xticklabels({'3D cases','Field piles','Field pairs'});ylabel('Count');ylim([0 5]);panel(3);
nexttile;hold on;plot([5 10 10 5 5],[2 2 6 6 2],':','Color',gray);scatter([7.5 8.75 8.75],[3 4 4],40,blue,'filled');scatter(geom.D_m,geom.L_installed_m./geom.D_m,40,orange,'filled');
xlabel('Diameter, {\it D} (m)');ylabel('{\it L/D}');xlim([1 11]);ylim([1 7]);legend({'Published geometry domain','Design cases','Field transfer'},'Location','southoutside','FontSize',7.5);panel(4);
saveall(f,'Fig04_evidence_coverage',4,'Geometry and independence of the validation evidence. (a) Diameter and embedment ratio of the two clay and two sand design cases, and four 2 m field piles. The sand cases share geometry but differ in relative density. (b) Corresponding wall thicknesses. (c) Four published 3D design conditions and four field specimens represent different evidence units; the field specimens form two simultaneously loaded, mutually reacting site pairs. (d) The published PISA geometric calibration domain and the field-scale extrapolation. Numerical design cases are distinct from the original calibration geometries but belong to the same model family.','Design geometry; pisa_2m_geometry.csv','Separates physical scale, calibration domain and independent evidence units.');

O=readtable(fullfile(root,'soil_profiles','offshore_profiles_0p05m.csv'));
CF=readtable(fullfile(root,'published_data','cowden_interpreted_soil_profiles.csv'));
DF=readtable(fullfile(root,'published_data','dunkirk_interpreted_G0_profile.csv'));
f=canvas(2,2,18,13);
nexttile;hold on;plot(O.su_kPa,O.z_m,'Color',blue);q=CF(string(CF.variable)=="su_kPa",:);plot(q.value,q.z_m,'o--','Color',green,'MarkerSize',3);set(gca,'YDir','reverse');ylim([0 20]);xlabel('Undrained strength, {\it s}_{u} (kPa)');ylabel('Depth, {\it z} (m)');legend({'Offshore reference','Cowden interpretation'},'Location','southwest');panel(1);
nexttile;hold on;plot(O.G0_MPa,O.z_m,'Color',blue);q=CF(string(CF.variable)=="G0_MPa",:);plot(q.value,q.z_m,'o--','Color',green,'MarkerSize',3);set(gca,'YDir','reverse');ylim([0 20]);xlabel('Shear modulus, {\it G}_{0} (MPa)');ylabel('Depth, {\it z} (m)');panel(2);
nexttile;hold on;zz=linspace(0,35,201)';for j=1:2,dr=[.55 .85];ss=pisa_dunkirk_soil(zz,dr(j));plot(ss.G0_Pa/1e6,zz,'Color',cols(j+2,:),'LineStyle',styles{j});end;plot(DF.value,DF.z_m,'o--','Color',gray,'MarkerSize',3);set(gca,'YDir','reverse');ylim([0 20]);xlabel('Shear modulus, {\it G}_{0} (MPa)');ylabel('Depth, {\it z} (m)');legend({'Offshore 55%','Offshore 85%','Dunkirk interpretation'},'Location','southwest');panel(3);
nexttile;hold on;errorbar(O.su_kPa(1:20:401),O.z_m(1:20:401),O.su_graphical_bound_kPa(1:20:401),'horizontal','.','Color',blue,'CapSize',2);set(gca,'YDir','reverse');ylim([0 20]);xlabel('Offshore {\it s}_{u} (kPa)');ylabel('Depth, {\it z} (m)');panel(4);
saveall(f,'Fig05_soil_inputs',4,'Distinct soil inputs used for the numerical and field comparisons. (a,b) Undrained strength and small-strain shear modulus for the representative offshore clay and the author-interpreted Cowden field profiles. (c) Representative offshore sand moduli at relative densities of 55% and 85%, compared with the interpreted onshore Dunkirk profile. (d) Conservative graphical reading ranges of the offshore strength path, including stroke width, depth sensitivity and axis registration. These ranges describe extraction resolution and are not soil variability or confidence intervals. The offshore reference profiles must not be substituted for measured or interpreted field profiles.','offshore_profiles_0p05m.csv; cowden_interpreted_soil_profiles.csv; dunkirk_interpreted_G0_profile.csv','Makes differences between hypothetical offshore and field soil inputs explicit.');

for group=1:2
 f=canvas(2,2,18,13);ci=(group-1)*2+(1:2);
 for ii=1:2
  c=cases(ci(ii));[~,gg]=design_config(c);
  for k=1:2
   if k==1,rg='ult';scale=1;xl='Ground displacement, {\it v}_{G} (m)';else,rg='sd';scale=1e3;xl='Ground displacement, {\it v}_{G} (mm)';end
   nexttile((k-1)*2+ii);hold on;
   [x,y]=published_curve(T,c,'3D',rg);plot(x*scale,y/1e6,'k-','LineWidth',1.1);
   [x,y]=published_curve(T,c,'published_1D',rg);plot(x*scale,y/1e6,'--','Color',gray);
   q=C(string(C.case_id)==c&string(C.model)=="Full",:);plot(q.head_w_m*scale,q.H_N/1e6,'Color',blue);
   if k==1,xlim([0 gg(1)/10]);else,xlim([0 gg(1)/10000*scale]);end
   xlabel(xl);ylabel('Lateral load, {\it H} (MN)');panel((k-1)*2+ii);
   if ii==1&&k==1,legend({'Published 3D','Published 1D','Four-reaction model'},'Location','southeast');end
  end
 end
 if group==1,id='Fig06_clay_design_curves';note='Columns show clay D1 (D=7.5 m, L=22.5 m, t=0.068 m, h=37.5 m) and D2 (D=8.75 m, L=35 m, t=0.091 m, h=87.5 m), respectively. Sources: Byrne et al. (2020), Fig. 17.';
 else,id='Fig07_sand_design_curves';note='Columns show sand D2t at relative densities 55% and 85%, respectively (D=8.75 m, L=35 m, t=0.150 m, h=87.5 m). Sources: Burd et al. (2020), Fig. 15.';end
 saveall(f,id,4,['Unfitted design-case load predictions against published numerical evidence. (a,b) Large-displacement paths to vG=D/10. (c,d) Small-displacement paths to vG=D/10000. ' note ' The published 3D and 1D paths are extracted from native PDF vectors. The present primary calculation uses E=200 GPa, nu=0.3, thin-wall properties and a Timoshenko shear factor of 0.5. Published curves are numerical predictions, not physical measurements.'],'predicted_curves.csv; published_design_curves.csv','Tests both initial stiffness and nonlinear capacity on non-calibration geometries.');
end

for group=1:2
 f=canvas(2,2,18,13);ci=(group-1)*2+(1:2);
 for ii=1:2
  c=cases(ci(ii));
  for k=1:2
   frac=[.5 1];frac=frac(k);nexttile((k-1)*2+ii);hold on;
   [z,M]=published_curve(T,c,'3D','m',frac);plot(M/1e6,z,'k-');
   [z,M]=published_curve(T,c,'published_1D','m',frac);plot(M/1e6,z,'--','Color',gray);
   q=P(string(P.case_id)==c&P.load_fraction==frac,:);plot(abs(q.M_theta_Nm)/1e6,q.z_m,'Color',blue);set(gca,'YDir','reverse');
   xlabel('Bending moment, {\it M} (MN m)');ylabel('Depth, {\it z} (m)');panel((k-1)*2+ii);
   if ii==1&&k==1,legend({'Published 3D','Published 1D','Four-reaction model'},'Location','northwest');end
  end
 end
 if group==1,id='Fig08_clay_moment_profiles';note='Columns show clay D1 and D2; original evidence is Byrne et al. (2020), Fig. 18.';else,id='Fig09_sand_moment_profiles';note='Columns show sand D2t at relative densities 55% and 85%; original evidence is Burd et al. (2020), Fig. 15.';end
 saveall(f,id,4,['Bending-moment distribution at matched externally defined loads. (a,b) H=0.5Href. (c,d) H=Href. Href is the published 3D load at ground displacement D/10, so every model is compared under the same force and eccentricity. ' note ' Moment magnitudes are shown; the implementation stores signed work-conjugate beam moments.'],'matched_load_moments.csv; published_design_curves.csv','Checks internal structural demand under common loading, beyond head-curve agreement.');
end

f=canvas(2,2,18,13);ee=E(string(E.reference)=="3D"&string(E.model)=="Full",:);
for k=1:2
 nexttile(k);hold on;rg=["sd","ult"];q=ee(string(ee.regime)==rg(k),:);
 errorbar(1:4,100*q.NMAE,100*(q.NMAE-q.NMAE_graphical_lower),100*(q.NMAE_graphical_upper-q.NMAE),'o','Color',blue,'MarkerFaceColor',blue,'CapSize',5);yline(10,':','Color',gray);xticks(1:4);xticklabels(labels);xlim([.5 4.5]);ylabel('Curve NMAE (%)');panel(k);
end
nexttile;hold on;for k=1:2,rg=["sd","ult"];q=ee(string(ee.regime)==rg(k),:);plot(1:4,100*(q.ratio_endpoint-1),'-o','Color',cols(k,:),'MarkerSize',4);end;yline(0,':k');xticks(1:4);xticklabels(labels);ylabel('Endpoint load error (%)');legend({'Small displacement','Large displacement'},'Location','southwest');panel(3);
nexttile;me=readtable(fullfile(d.run,'moment_metrics.csv'));vals=zeros(4,2);for j=1:4,vals(j,:)=100*me.moment_NMAE(string(me.case_id)==cases(j))';end;bb=bar(vals);bb(1).FaceColor=orange;bb(2).FaceColor=green;xticks(1:4);xticklabels(labels);ylabel('Moment-profile NMAE (%)');legend({'0.5{\it H}_{ref}','{\it H}_{ref}'},'Location','northeast');panel(4);
saveall(f,'Fig10_prediction_errors',4,'Prediction error resolved by case and response regime. (a,b) Load-curve normalized mean absolute error over small- and large-displacement domains, respectively. Bars propagate native-curve stroke width and horizontal coordinate resolution; they are not confidence intervals. The 10% line is a prespecified research screening value, not a design-code limit. (c) Signed load errors at D/10000 and D/10. (d) Bending-moment profile errors at 0.5 and 1.0 times the published 3D reference load at ground displacement D0/10. Each design condition is reported separately; the number of digitized coordinates does not increase independent evidence.','curve_metrics.csv; moment_metrics.csv','Distinguishes stiffness, capacity and internal-force accuracy and graphical resolution.');

f=canvas(2,2,18,13);bf=B(string(B.model)=="Full",:);mom=[bf.p_resultant_moment_Nm,bf.distributed_moment_Nm,bf.base_H_moment_Nm,bf.base_M_Nm];mom=mom./sum(mom,2);
reactionTickLabels={'Clay\newlineD1','Clay\newlineD2','Sand\newline55%','Sand\newline85%'};
nexttile;bb=bar(mom,'stacked');for j=1:4,bb(j).FaceColor=cols(j,:);end;xticks(1:4);xticklabels(reactionTickLabels);ylabel({'Signed moment','equilibrium fraction'});legend({'{\it p}','{\it m}','{\it H}_{B}','{\it M}_{B}'},'Location','southoutside','Orientation','horizontal');panel(1);
nexttile;forces=[bf.distributed_H_N bf.base_H_N];forces=forces./sum(forces,2);bb=bar(forces,'stacked');bb(1).FaceColor=blue;bb(2).FaceColor=green;xticks(1:4);xticklabels(reactionTickLabels);ylabel({'Signed force','equilibrium fraction'});legend({'Distributed force','Base shear'},'Location','southoutside','Orientation','horizontal');panel(2);
UP=readtable(fullfile(d.run,'ultimate_profiles.csv'));
for k=1:2,nexttile(k+2);hold on;for j=1:4,q=UP(string(UP.case_id)==cases(j)&string(UP.model)=="Full",:);[~,gg]=design_config(cases(j));if k==1,xx=q.p_N_per_m/1e6;else,xx=q.m_N/1e6;end;plot(xx,q.z_m/gg(3),'Color',cols(j,:),'LineStyle',styles{j});end;set(gca,'YDir','reverse');ylabel('Normalized depth, {\it z/L}');if k==1,xlabel('Distributed force, {\it p} (MN/m)');else,xlabel('Distributed moment, {\it m} (MN)');end;panel(k+2);end
legend(labels,'Location','southoutside','NumColumns',2);
saveall(f,'Fig11_reaction_channels',4,'Four-reaction load transfer at each model ground displacement D/10. (a) Signed shares of moment equilibrium from the moment of the distributed force, distributed moment, base-shear lever arm and base moment. (b) Signed horizontal-force shares. Opposing reactions may produce negative or greater-than-unity shares; absolute-value renormalization is not used. (c,d) Distributed lateral force and distributed moment along normalized depth. Different cases have different loads at D/10. These are internal model reaction partitions and do not constitute independently observed component forces.','reaction_balance.csv; ultimate_profiles.csv','Shows where force and moment equilibrium are carried and avoids confusing m units with p.');

mods=["Full","p_only","No_m","No_HB","No_MB"];
f=canvas(2,2,18,13);
for j=1:4,nexttile;hold on;values=zeros(5,2);for k=1:2,rg=["sd","ult"];for im=1:5,q=E(string(E.case_id)==cases(j)&string(E.model)==mods(im)&string(E.reference)=="3D"&string(E.regime)==rg(k),:);values(im,k)=q.ratio_endpoint;end;end;bb=bar(values);bb(1).FaceColor=blue;bb(2).FaceColor=orange;yline(1,':k');xticks(1:5);xticklabels({'Full','Only {\it p}','No {\it m}','No {\it H}_{B}','No {\it M}_{B}'});ylabel('Predicted / 3D endpoint load');panel(j);if j==1,legend({'Small displacement','Large displacement'},'Location','southoutside','Orientation','horizontal');end;end
saveall(f,'Fig12_component_ablation',4,'Dependence on the four calibrated reaction channels. Panels (a–d) show clay D1, clay D2, sand 55% and sand 85%, respectively. Endpoint load ratios are reported at D/10000 and D/10 for the full model, the p-only model and removal of one additional reaction channel at a time. All surviving coefficients remain unchanged. Because second-stage calibration can compensate among components, these ablations measure internal functional dependence of this calibrated representation and are not direct estimates of experimental causal contributions.','curve_metrics.csv','Tests whether omitting rotational and toe reactions can preserve both stiffness and capacity.');

f=canvas(2,2,18,13);
for j=1:4,nexttile;hold on;mm=["Full","EB_limit","Exact_annulus"];for k=1:3,q=C(string(C.case_id)==cases(j)&string(C.model)==mm(k),:);plot(q.head_w_m*1e3,q.H_N/1e6,'Color',cols(k,:),'LineStyle',styles{k});end;[~,gg]=design_config(cases(j));xlim([0 gg(1)/10000*1e3]);xlabel('Ground displacement, {\it v}_{G} (mm)');ylabel('Lateral load, {\it H} (MN)');panel(j);if j==1,legend({'Thin wall, Timoshenko','Euler–Bernoulli limit','Exact annulus, Timoshenko'},'Location','southeast');end;end
saveall(f,'Fig13_beam_assumptions',4,'Structural representation sensitivity in the small-displacement regime. Panels (a–d) follow the four design cases in the established order. The primary thin-wall Timoshenko representation is compared with a numerical Euler–Bernoulli limit (shear factor 5000) and exact annular section properties, retaining the same soil reactions and steel modulus. These changes are declared alternatives and are not selected by the best fit to the published curves. Their effects help interpret small differences among reported implementations.','predicted_curves.csv','Separates shear deformation and section conventions from soil-model discrepancies.');

for view=1:2
 f=canvas(2,2,18,13);
 for j=1:4
  nexttile;hold on;for k=1:3,frac=[.1 .5 1];frac=frac(k);u=G(string(G.case_id)==cases(j)&G.load_fraction==frac&string(G.treatment)=="Updated",:);fr=G(string(G.case_id)==cases(j)&G.load_fraction==frac&string(G.treatment)=="Frozen",:);
   if view==1
    base=u.w_m(u.D_ratio==1);plot(u.D_ratio,u.w_m/base,'Color',blue,'LineStyle',styles{k},'Marker','o','MarkerSize',3);plot(fr.D_ratio,fr.w_m/base,'Color',orange,'LineStyle',styles{k},'Marker','s','MarkerSize',3);
    if ismember('actual_geometry_domain',u.Properties.VariableNames),bad=~u.actual_geometry_domain|~u.w_within_0p1D;plot(u.D_ratio(bad),u.w_m(bad)/base,'kx','MarkerSize',5,'HandleVisibility','off');bad=~fr.actual_geometry_domain|~fr.w_within_0p1D;plot(fr.D_ratio(bad),fr.w_m(bad)/base,'kx','MarkerSize',5,'HandleVisibility','off');end
   else
    bias=100*(fr.w_m./u.w_m-1);plot(u.D_ratio,bias,'Color',cols(k,:),'LineStyle',styles{k},'Marker','o','MarkerSize',3);
    if ismember('actual_geometry_domain',u.Properties.VariableNames),bad=~u.actual_geometry_domain|~u.w_within_0p1D|~fr.w_within_0p1D;plot(u.D_ratio(bad),bias(bad),'kx','MarkerSize',5,'HandleVisibility','off');end
   end
  end
  xlabel('Actual / nominal diameter, {\it D/D}_{0}');if view==1,ylabel('Displacement / nominal value');else,ylabel('Frozen-diameter displacement bias (%)');yline(0,':','Color',gray);end;panel(j);
  if view==1,yl=ylim;ylim([yl(1),yl(2)+.16*diff(yl)]);xlim([.89 1.12]);end
  if j==1,if view==1,lg=legend({'Updated, 0.1','Frozen, 0.1','Updated, 0.5','Frozen, 0.5','Updated, 1.0','Frozen, 1.0'},'NumColumns',3,'FontSize',8);lg.Layout.Tile='south';else,legend({'0.1{\it H}_{ref}','0.5{\it H}_{ref}','{\it H}_{ref}'},'Location','northwest');end;end
 end
 if view==1,id='Fig14_geometry_counterfactuals';caption='Matched-load geometric counterfactuals. Panels (a–d) show clay D1, clay D2, sand 55% and sand 85%. Actual pile diameter changes by up to 10%, while wall thickness, embedment, eccentricity and material inputs remain fixed. In both treatments the beam uses the actual diameter; only the soil-reaction diameter is either updated or fixed at its nominal value. Line styles identify 0.1, 0.5 and 1.0 times the original published 3D reference load at ground displacement D0/10. Each path is normalized by its own nominal response. High-load perturbations can exceed the original displacement calibration range and remain flagged in the data.';purpose='Separates structural and soil geometric updates under identical external actions.';
 else,id='Fig15_frozen_diameter_bias';caption='Signed displacement bias caused by freezing the soil-reaction diameter. Panels (a–d) follow the four design cases. Bias is 100 times the frozen-to-updated displacement ratio minus one at common actual geometry and load. The nominal-diameter null control is exactly zero. Load fractions refer to the original 3D reference load at ground displacement D0/10; they are deterministic interventions, not event probabilities. These error paths are conditional on the reduced model and do not independently validate its derivative against a perturbed 3D or field experiment.';purpose='Shows the sign, magnitude and load dependence of geometry-freezing error.';end
 caption=[caption ' Black crosses flag actual geometry or displacement outside the published calibration bounds.'];saveall(f,id,4,caption,'geometry_counterfactuals.csv',purpose);
end

f=canvas(2,2,18,13);
for j=1:4,nexttile;hold on;q=A(string(A.case_id)==cases(j)&A.relative_step==1e-4,:);bb=bar(q.load_fraction,-[q.beam_derivative q.soil_derivative],'stacked');bb(1).FaceColor=blue;bb(2).FaceColor=orange;plot(q.load_fraction,-q.central_derivative,'ks','MarkerSize',5,'LineStyle','none');xticks([.1 .5 1]);xlabel('Load ratio, {\it H}/{\it H}_{ref}');ylabel({'Negative displacement','derivative (m/m)'});panel(j);if j==1,legend({'Beam term','Soil term','Central difference'},'Location','northwest');end;end
saveall(f,'Fig16_adjoint_decomposition',4,'Work-conjugate decomposition of the diameter derivative. Panels (a–d) follow the four design cases. Bars separate the beam and soil contributions to minus dv_G/dD at fixed load, thickness, embedment and eccentricity; symbols are central finite-difference derivatives. The transpose system J^T λ = e_v is solved for the full, possibly nonsymmetric, tangent. A relative diameter step of 10^-4 is used for this displayed cross-check; a second step is retained in the verification data. No experimental causal interpretation is assigned to the algebraic partition.','adjoint_derivatives.csv','Identifies how soil geometry becomes increasingly important as resistance is mobilized.');

f=canvas(2,4,18,16);zd=table;
for j=1:4
 q=Z(string(Z.case_id)==cases(j)&Z.load_fraction==1,:);[cfg,gg]=design_config(cases(j));aa=pisa_assemble(cfg,gg(1),gg(2),gg(3));dp=q.p_derivative_contribution./aa.weights;dm=q.m_derivative_contribution./aa.weights;
 nexttile(j);hold on;plot(dp,q.z_m/gg(3),'Color',blue);plot(dm,q.z_m/gg(3),'Color',orange);set(gca,'YDir','reverse');ax=gca;ax.XAxis.Exponent=0;xtickformat('%.3f');xlabel({'Local sensitivity','density (m^{-1})'});ylabel('Normalized depth, {\it z/L}');panel(j);
 nexttile(j+4);hold on;plot(cumsum(q.p_derivative_contribution),q.z_m/gg(3),'Color',blue);plot(cumsum(q.m_derivative_contribution),q.z_m/gg(3),'Color',orange);plot(cumsum(q.p_derivative_contribution+q.m_derivative_contribution),q.z_m/gg(3),'--','Color',green);set(gca,'YDir','reverse');xlabel({'Cumulative sensitivity','(m/m)'});ylabel('Normalized depth, {\it z/L}');panel(j+4);
 if j==1,lg=legend({'{\it p} term','{\it m} term','Distributed sum'},'NumColumns',3);lg.Layout.Tile='south';end
 zq=table(repmat(cases(j),height(q),1),q.z_m,dp,dm,cumsum(q.p_derivative_contribution),cumsum(q.m_derivative_contribution),'VariableNames',{'case_id','z_m','p_density_per_m','m_density_per_m','p_cumulative','m_cumulative'});zd=[zd;zq];
end
writetable(zd,fullfile(out,'Fig17_plot_data.csv'));
for ax=reshape(findall(f,'Type','axes'),1,[]),ax.YLim=[-.07 1.02];ax.YTick=0:.2:1;ax.YDir='reverse';end
for ax=reshape(findall(f,'Type','axes'),1,[]),ax.FontSize=10.5;ax.XLabel.FontSize=10.5;ax.YLabel.FontSize=10.5;end
for lg=reshape(findall(f,'Type','legend'),1,[]),lg.FontSize=10;end
saveall(f,'Fig17_depth_resolved_sensitivity',8,'Depth-resolved soil sensitivity at the published 3D reference load at ground displacement D0/10. Columns follow clay D1, clay D2, sand 55% and sand 85%. (a–d) Local adjoint-weighted diameter-derivative densities from distributed force and moment, obtained by dividing each quadrature contribution by its integration weight. (e–h) Separate and combined cumulative distributed contributions. Each row has its own dimensional scale so local contributions remain visible. Signs are preserved. Base contributions are excluded from the cumulative curves and are reported separately in the complete four-channel accounting. No normalization as feature importance or experimental causal interpretation is used.','depth_derivative_contributions.csv; Fig17_plot_data.csv','Separates local influence density from accumulated response sensitivity without hiding one by a shared scale.');

f=canvas(1,2,18,7.8);
nexttile;hold on;for j=1:4,q=A(string(A.case_id)==cases(j)&A.relative_step==1e-4,:);semilogy(q.load_fraction,max(q.relative_difference,1e-15),'-o','Color',cols(j,:),'MarkerSize',4);end;set(gca,'YScale','log');xlabel('Load ratio, {\it H}/{\it H}_{ref}');ylabel('Relative derivative discrepancy');legend(labels,'Location','northwest','FontSize',7.5);panel(1);
nexttile;hold on;for j=1:4,q=A(string(A.case_id)==cases(j)&A.relative_step==1e-4,:);plot(q.load_fraction,100*(q.incorrect_no_transpose_derivative./q.adjoint_derivative-1),'-o','Color',cols(j,:),'MarkerSize',4);end;yline(0,':','Color',gray);xlabel('Load ratio, {\it H}/{\it H}_{ref}');ylabel('Error without adjoint transpose (%)');panel(2);
saveall(f,'Fig18_adjoint_consistency',2,'Numerical and algebraic consistency of the geometry derivative. (a) Relative discrepancy between the transpose-adjoint derivative and a central finite difference with relative step 10^-4. (b) Error produced by incorrectly solving J λ = e_v instead of J^T λ = e_v. The clay implementation has a symmetric tangent, whereas the sand distributed moment depends on lateral reaction and introduces a cross derivative dm/dv. Small numerical agreement in panel a verifies implementation of the model derivative; it does not verify a physical mechanism independently.','adjoint_derivatives.csv','Demonstrates why work conjugacy and a nonsymmetric adjoint matter for sand.');

f=canvas(2,2,18,13);scenarios=["Nominal","G0_0p8","G0_1p2","Strength_low","Strength_high"];
for j=1:4,nexttile;hold on;cc=[.15 .15 .15;blue;blue;orange;orange];st={'-','--','-','--','-'};for k=1:5,q=S(string(S.case_id)==cases(j)&string(S.scenario)==scenarios(k),:);plot(q.head_w_m,q.H_N/1e6,'Color',cc(k,:),'LineStyle',st{k});end;xlabel('Ground displacement, {\it v}_{G} (m)');ylabel('Lateral load, {\it H} (MN)');panel(j);if j==1,legend({'Nominal','0.8{\it G}_{0}','1.2{\it G}_{0}','Lower strength','Higher strength'},'Location','southeast','FontSize',7.5);elseif j==3,legend({'Nominal','0.8{\it G}_{0}','1.2{\it G}_{0}','Lower relative density','Higher relative density'},'Location','southeast','FontSize',7.5);end;end
saveall(f,'Fig19_soil_sensitivity',4,'Prespecified soil-property scenarios without response fitting. Panels (a–d) follow the four design cases. Small-strain modulus is multiplied by 0.8 or 1.2. Strength scenarios multiply the clay undrained strength by 0.8 or 1.2, or vary sand relative density by minus or plus 0.05 with the corresponding published modulus relation. Other inputs remain fixed. These deterministic ranges are sensitivity scenarios, not estimates of site uncertainty or confidence limits; no best-fitting scenario is substituted for the nominal prediction.','soil_sensitivity_curves.csv','Distinguishes stiffness and capacity responses to independently specified parameter changes.');

V=readtable(fullfile(root,'pisa_model','results','published_sand_table7_verification.csv'));
f=canvas(1,2,18,7.8);
for k=1:2,nexttile;hold on;for j=1:2,ids=["C1","C4"];q=V(string(V.case_id)==ids(j),:);[nn,ord]=sort(q.elements);if k==1,vals=abs(q.Hult_MN-q.Hult_published_MN);else,vals=abs(q.Hsd_MN-q.Hsd_published_MN);end;loglog(nn,max(vals(ord),1e-9),'-o','Color',cols(j,:),'MarkerSize',4);end;set(gca,'XScale','log','YScale','log');yline(.00005,':k');xlabel('Number of embedded elements');ylabel('Published-table difference (MN)');panel(k);if k==1,lg=legend({'C1','C4','Half last reported digit'},'Location','southoutside','Orientation','horizontal');lg.Layout.Tile='south';end;end
saveall(f,'Fig20_published_table_reproduction',2,'Direct reproduction of the published sand mesh study. (a,b) Absolute differences in ultimate and small-displacement loads, respectively, from Burd et al. (2020), Table 7, using each original mesh. The dotted line is half the last reported 0.0001 MN digit. Thirteen case-and-mesh combinations are compared without response calibration. The original published convergence sequence is being reproduced; closeness to its rounded numbers is numerical implementation evidence, not a new physical validation. Explicit fine-mesh solver tolerances and all failed preliminary checks are retained.','published_sand_table7_verification.csv','Anchors implementation accuracy to exact published tabulated values, independent of plot extraction.');

writetable(cell2table(reg,'VariableNames',{'id','panels','caption','data_sources','evidence_value','program'}),fullfile(out,'figure_captions_evidence.csv'));
fprintf('Exported %d vector group figures\n',size(reg,1));
 function f=canvas(nr,nc,w,h),f=figure('Visible','off','Units','centimeters','Position',[1 1 w h],'Color','w','Renderer','painters');tiledlayout(nr,nc,'TileSpacing','compact','Padding','compact');end
 function panel(k),title('');text(.015,.985,['(' char('a'+k-1) ')'],'Units','normalized','VerticalAlignment','top','FontSize',10,'FontWeight','bold','Tag','PanelLabel');end
 function saveall(f,id,np,caption,data,purpose)
  if ~isempty(selected)&&~any(startsWith(string(id),selected)),reg(end+1,:)={id,np,caption,data,purpose,'code/plot_validation_figures.m'};close(f);return;end
  for ax=reshape(findall(f,'Type','axes'),1,[]),ax.FontSize=9;ax.XLabel.FontSize=10;ax.YLabel.FontSize=10;end
  for lg=reshape(findall(f,'Type','legend'),1,[]),lg.FontSize=8.5;end
  for tx=reshape(findall(f,'Tag','PanelLabel'),1,[]),tx.FontSize=10.5;end
  drawnow;axs=findall(f,'Type','axes');for ax=axs',ax.Title.String='';ax.Box='off';end
  if ismember(id,{'Fig05_soil_inputs','Fig10_prediction_errors'}),f.PaperUnits='centimeters';f.Units='centimeters';fp=f.Position;f.PaperPosition=[0 0 fp(3:4)];f.PaperSize=fp(3:4);f.PaperPositionMode='manual';annotation(f,'rectangle',[0.002 0.002 0.996 0.996],'Color',[1 1 1],'LineWidth',0.1,'Tag','PublicationCanvasBoundary');drawnow;end
  savefig(f,fullfile(out,[id '.fig']));exportgraphics(f,fullfile(out,[id '.pdf']),'ContentType','vector','BackgroundColor','white');print(f,fullfile(out,[id '.svg']),'-dsvg');print(f,fullfile(out,[id '.emf']),'-dmeta');exportgraphics(f,fullfile(out,[id '_QA.png']),'Resolution',180);
  reg(end+1,:)={id,np,caption,data,purpose,'code/plot_validation_figures.m'};close(f);
 end
end

