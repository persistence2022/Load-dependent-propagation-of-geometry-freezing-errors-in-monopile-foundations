function plot_additional_checks
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'code'),fullfile(root,'pisa_model'));out=fullfile(root,'figures');
[~,mi]=resolve_current_run(fullfile(root,'current_mechanism.json'),root);[~,di]=resolve_current_run(fullfile(root,'current_design_formal.json'),root);[~,fi]=resolve_current_run(fullfile(root,'field_validation','current_formal.json'),root);[~,bi]=resolve_current_run(fullfile(root,'design_baselines','current_formal.json'),root);
A=readtable(fullfile(mi.run,'adjoint_derivatives.csv'));B=readtable(fullfile(root,'pisa_model','results','independent_bvp_matlab_comparison.csv'));F=readtable(fullfile(fi.run,'numerical_checks','mesh_convergence.csv'));
C=readtable(fullfile(bi.run,'predicted_curves.csv'));E=readtable(fullfile(bi.run,'curve_metrics.csv'));PC=readtable(fullfile(di.run,'predicted_curves.csv'));PE=readtable(fullfile(di.run,'curve_metrics.csv'));
T=readtable(fullfile(root,'sources','design_benchmarks','published_design_curves.csv'));cases=["Clay_D1","Clay_D2","Sand_D2t_Dr55","Sand_D2t_Dr85"];
blue=[0 114 178]/255;orange=[213 94 0]/255;green=[0 158 115]/255;purple=[204 121 167]/255;colors=[blue;orange;green;purple];reg={};
set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman','defaultAxesFontSize',9,'defaultTextFontSize',9,'defaultAxesLineWidth',.7,'defaultLineLineWidth',1.25,'defaultAxesTickDir','out','defaultLegendBox','off','defaultLegendFontSize',8);
f=canvas(2,2,18,13);
for j=1:4,nexttile;hold on;q=A(string(A.case_id)==cases(j)&A.relative_step==1e-4,:);v=[q.p_derivative,q.m_derivative,q.HB_derivative,q.MB_derivative];v=v./q.soil_derivative;bb=bar(q.load_fraction,v,'stacked');for k=1:4,bb(k).FaceColor=colors(k,:);end;xticks([.1 .5 1]);xlabel('Load ratio, {\it H}/{\it H}_{ref}');ylabel('Signed soil-derivative fraction');ylim([0 1.2]);panel(j);if j==1,legend({'{\it p}','{\it m}','{\it H}_{B}','{\it M}_{B}'},'Location','southoutside','Orientation','horizontal');end;end
saveall(f,'Fig27_soil_derivative_channels',4,'Complete four-channel accounting of the soil-diameter derivative. Panels (a–d) show clay D1, clay D2, sand 55% and sand 85%. Each signed contribution is divided by the total soil derivative at the same load. The base-shear and base-moment terms are included explicitly, complementing the distributed depth contributions. Closure against the full residual derivative is checked numerically. The partition belongs to the stated calibrated model and is not a direct experimental separation of soil mechanisms.','adjoint_derivatives.csv','Closes the distributed and toe sensitivity budget rather than omitting end reactions.');

f=canvas(2,2,18,13);
for k=1:2
 nexttile(k);hold on;for j=1:2,L=[20 60];q=B(B.L_m==L(j),:);nn=unique(q.elements);val=zeros(size(nn));for ii=1:numel(nn),if k==1,val(ii)=max(q.H_relative_error(q.elements==nn(ii)));else,val(ii)=max(q.theta_relative_error(q.elements==nn(ii)));end;end;loglog(nn,100*val,'-o','Color',colors(j,:),'MarkerSize',4);end;set(gca,'XScale','log','YScale','log');yline(.1,':k');xlabel('Number of embedded elements');if k==1,ylabel('Maximum load error (%)');else,ylabel('Maximum rotation error (%)');end;panel(k);if k==1,legend({'Sand C1','Sand C4','0.1% screening'},'Location','northeast');end
end
ids=unique(strcat(string(F.pile_id),"/",string(F.scenario_id)),'stable');
for k=1:2,nexttile(k+2);hold on;for j=1:numel(ids),mask=strcat(string(F.pile_id),"/",string(F.scenario_id))==ids(j);q=F(mask,:);if k==1,v=abs(q.Kinit_rel_to120);else,v=abs(q.H01D_rel_to120);end;plot(q.elements(q.elements<120),100*v(q.elements<120),'-o','Color',colors(j,:),'MarkerSize',4);end;set(gca,'YScale','log');xlabel('Number of embedded elements');if k==1,ylabel('Stiffness change to 120 elements (%)');else,ylabel('Load change to 120 elements (%)');end;panel(k+2);end
legend({'DL1, PISA','DL1, API 36°','CL1, PISA','CL1, Matlock'},'Interpreter','none','Location','southoutside','NumColumns',2,'FontSize',8);
axs=findall(f,'Type','axes');for ax=axs',yl=ylim(ax);ylim(ax,[yl(1)*.6 yl(2)*4]);if strcmp(ax.XScale,'log'),xlim(ax,[16 235]);xticks(ax,[20 40 80 200]);else,xlim(ax,[35 83]);xticks(ax,[40 80]);end;end
saveall(f,'Fig28_independent_numerical_checks',4,'Independent solution and mesh checks. (a,b) Maximum load and head-rotation discrepancies over five displacement targets per geometry between MATLAB finite elements and a separately implemented adaptive collocation solution of the four first-order beam equations. Both solve the same constitutive model. The retained 20-element C4 load discrepancy exceeds the prespecified 0.1% screening line; default 40-element and finer cases pass. (c,d) Stiffness and D/10-load differences from the 120-element field reference for four representative model/specimen configurations. Zero reference differences at 120 elements are omitted on the logarithmic axes. Mesh checks verify numerical implementation, not physical truth.','independent_bvp_matlab_comparison.csv; field numerical_checks/mesh_convergence.csv','Separates independent discretization agreement, mesh convergence and physical validation.');

f=canvas(2,2,18,13);
for j=1:4
 nexttile;hold on;[x,y]=published_curve(T,cases(j),'3D','ult');plot(x,y/1e6,'k-');
 q=PC(string(PC.case_id)==cases(j)&string(PC.model)=="Full",:);plot(q.head_w_m,q.H_N/1e6,'Color',blue);
 mods=unique(string(C.model(string(C.case_id)==cases(j))),'stable');styles={'--','-.',':'};
 for k=1:3,q=C(string(C.case_id)==cases(j)&string(C.model)==mods(k),:);plot(q.head_w_m,q.H_N/1e6,'Color',orange,'LineStyle',styles{k});end
 [~,g]=design_config(cases(j));xlim([0 g(1)/10]);xlabel('Ground displacement, {\it v}_{G} (m)');ylabel('Lateral load, {\it H} (MN)');panel(j);
 if j==1,legend({'Published 3D','Four reactions','Baseline low','Baseline middle','Baseline high'},'Location','southeast','FontSize',7.5);end
end
saveall(f,'Fig29_large_pile_baseline_curves',4,'Conventional p-y scenario predictions on the same large-pile design inputs. Panels (a–d) show clay D1, clay D2, sand 55% and sand 85%. The three baseline levels correspond to Matlock epsilon50=0.004, 0.005, 0.007 at J=0.5 for clay, and API peak friction angles 32, 36, 40 degrees for sand. Geometry, soil profiles, E=200 GPa, thin-wall properties and Timoshenko beam theory match the primary model. Baseline constitutive inputs are declared assumptions, not values inferred from reference curves. For sand, equal assumed friction angles yield equal API curves in the two density cases because a site-specific density-to-peak-angle relation is unavailable; this does not establish that API cannot represent density effects.','design_baselines/predicted_curves.csv; design predicted_curves.csv; published_design_curves.csv','Tests the original model family externally while exposing additional-input assumptions.');

f=canvas(2,2,18,13);
for j=1:4
 nexttile;hold on;plot([-10 10 10 -10 -10],[-10 -10 10 10 -10],':','Color',[.6 .6 .6]);plot(0,0,'k+','MarkerSize',7);
 q=PE(string(PE.case_id)==cases(j)&string(PE.model)=="Full"&string(PE.reference)=="3D",:);
 plot(100*(q.ratio_endpoint(string(q.regime)=="sd")-1),100*(q.ratio_endpoint(string(q.regime)=="ult")-1),'o','Color',blue,'MarkerFaceColor',blue,'MarkerSize',5);
 mods=unique(string(E.model(string(E.case_id)==cases(j))),'stable');v=zeros(3,2);
 for k=1:3,q=E(string(E.case_id)==cases(j)&string(E.model)==mods(k),:);v(k,:)=[100*(q.ratio_endpoint(string(q.regime)=="sd")-1),100*(q.ratio_endpoint(string(q.regime)=="ult")-1)];end
 plot(v(:,1),v(:,2),'s--','Color',orange,'MarkerSize',5);
 xlabel('Small-displacement endpoint error (%)');ylabel('Large-displacement endpoint error (%)');panel(j);if j==1,legend({'10% screening window','Published 3D','Four reactions','Baseline scenarios'},'Location','northeast','FontSize',7.5);end
end
saveall(f,'Fig30_baseline_regime_tradeoffs',4,'Joint small- and large-displacement errors in the external numerical comparison. Panels (a–d) follow the four design cases. The blue marker is the fixed primary four-reaction prediction; the orange path contains all three predeclared baseline parameter levels. Axes use signed errors at D/10000 and D/10, respectively. The origin denotes exact agreement with the published 3D reference and the square is the prespecified 10% research screening window. A near-capacity match need not imply a stiffness match. Parameter selection by minimum error is not performed, and the numerical reference is not described as field truth.','design_baselines/curve_metrics.csv; design curve_metrics.csv','Tests whether apparent capacity agreement also predicts the initial response.');
writetable(cell2table(reg,'VariableNames',{'id','panels','caption','data_sources','evidence_value','program'}),fullfile(out,'additional_figure_captions_evidence.csv'));
 function f=canvas(nr,nc,w,h),f=figure('Visible','off','Units','centimeters','Position',[1 1 w h],'Color','w','Renderer','painters');tiledlayout(nr,nc,'TileSpacing','compact','Padding','compact');end
 function panel(k),title('');text(.015,.985,['(' char('a'+k-1) ')'],'Units','normalized','VerticalAlignment','top','FontSize',10,'FontWeight','bold','Tag','PanelLabel');end
 function saveall(f,id,np,caption,data,purpose)
  for ax=reshape(findall(f,'Type','axes'),1,[]),ax.FontSize=9;ax.XLabel.FontSize=10;ax.YLabel.FontSize=10;end
  for lg=reshape(findall(f,'Type','legend'),1,[]),lg.FontSize=8.5;end
  for tx=reshape(findall(f,'Tag','PanelLabel'),1,[]),tx.FontSize=10.5;end
  drawnow;for ax=reshape(findall(f,'Type','axes'),1,[]),ax.Box='off';end
  savefig(f,fullfile(out,[id '.fig']));exportgraphics(f,fullfile(out,[id '.pdf']),'ContentType','vector');print(f,fullfile(out,[id '.svg']),'-dsvg');print(f,fullfile(out,[id '.emf']),'-dmeta');exportgraphics(f,fullfile(out,[id '_QA.png']),'Resolution',180);reg(end+1,:)={id,np,caption,data,purpose,'code/plot_additional_checks.m'};close(f);
 end
end
