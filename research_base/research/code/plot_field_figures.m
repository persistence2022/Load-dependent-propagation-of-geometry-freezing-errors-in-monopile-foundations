function plot_field_figures
root=fileparts(fileparts(mfilename('fullpath')));out=fullfile(root,'figures');
addpath(fullfile(root,'code'));[run,~]=resolve_current_run(fullfile(root,'field_validation','current_formal.json'),root);
C=readtable(fullfile(run,'predicted_curves.csv'));E=readtable(fullfile(run,'field_metrics.csv'));P=readtable(fullfile(run,'ultimate_profiles.csv'));
T=readtable(fullfile(root,'published_data','pisa_2m_published_paths.csv'));O=readtable(fullfile(root,'published_data','pisa_2m_response_metrics.csv'));
blue=[0 114 178]/255;orange=[213 94 0]/255;green=[0 158 115]/255;purple=[204 121 167]/255;gray=[.65 .65 .65];cols=[blue;orange;green;purple];
set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman','defaultAxesFontSize',9,'defaultTextFontSize',9,'defaultAxesLineWidth',.7,'defaultLineLineWidth',1.25,'defaultAxesTickDir','out','defaultLegendFontSize',8,'defaultLegendBox','off');reg={};
for group=1:2
 if group==1,piles=["DL1","DL2"];site='Dunkirk';id='Fig21_dunkirk_field_curves';else,piles=["CL1","CL2"];site='Cowden';id='Fig22_cowden_field_curves';end
 f=canvas(2,2,18,13);
 for j=1:2
  for k=1:2
   nexttile((k-1)*2+j);hold on;
   if k==1,pan="a";xmax=.2;scale=1;xl='Ground displacement, {\it v}_{G} (m)';else,pan="b";xmax=.004;scale=1e3;xl='Ground displacement, {\it v}_{G} (mm)';end
   q=T(string(T.pile_id)==piles(j)&string(T.panel)==pan&T.inside_axes==1,:);
   groups=unique(q.path_index,'stable');for jj=1:numel(groups),s=q(q.path_index==groups(jj),:);ph=plot(s.v_ground_m*scale,s.H_N/1e6,'Color',gray,'LineWidth',.65);end
   a=C(string(C.pile_id)==piles(j)&string(C.scenario_id)=="PISA_full",:);lh1=plot(a.head_w_m*scale,a.H_N/1e6,'Color',blue);
   a=C(string(C.pile_id)==piles(j)&string(C.scenario_id)=="PISA_p_only",:);lh2=plot(a.head_w_m*scale,a.H_N/1e6,'--','Color',orange);
   obs=O(string(O.pile_id)==piles(j),:);if k==1,xx=.2;yy=obs.H_0p1D_kN*1e3;else,xx=.002;yy=obs.H_at_initial_endpoint_N_derived;end
   oh=plot(xx*scale,yy/1e6,'kd','MarkerFaceColor','w','MarkerSize',5,'LineWidth',1.1);xlim([0 xmax*scale*1.025]);publishedMax=max(q.H_N(q.v_ground_m>=0&q.v_ground_m<=xmax),[],'omitnan')/1e6;
   ylim([0 max([yy/1e6,interp1(a.head_w_m,a.H_N,xmax)/1e6,publishedMax])*1.12]);
   xlabel(xl);ylabel('Lateral load, {\it H} (MN)');panel((k-1)*2+j);
   if j==1&&k==1,legend([ph lh1 lh2 oh],{'Published test history','Four-reaction model','{\it p}-only model','Published Table 3'},'Location','southeast');end
  end
 end
 cap=['Unfitted comparison with the ' site ' 2 m field piles. Columns show ' char(piles(1)) ' and ' char(piles(2)) '. (a,b) Ground load-displacement response to D/10. (c,d) An expanded initial response with the author-interpreted Table 3 point at D/1000. Gray paths retain published loading, unloading and hold-related history where visible; they are not fitted or converted into independent samples. Diamonds are the primary author-processed response metrics. Monotonic model curves use the site-specific interpreted soil profiles and explicit extrapolation beyond the PISA diameter calibration domain. The two piles were simultaneously loaded against one another and constitute one site-pair unit.'];
 saveall(f,id,4,cap,'predicted_curves.csv; pisa_2m_published_paths.csv; pisa_2m_response_metrics.csv','Displays the two primary field metrics together with the correlated published loading history.');
end

piles=["DL1","DL2","CL1","CL2"];f=canvas(1,2,18,8.3);
for k=1:2,nexttile;hold on;v=zeros(4,2);for j=1:4,for mm=1:2,models=["PISA_full","PISA_p_only"];q=E(string(E.pile_id)==piles(j)&string(E.scenario_id)==models(mm),:);if k==1,v(j,mm)=q.Kinit_relative_error;else,v(j,mm)=q.H01D_relative_error;end;end;end;bb=bar(100*v);bb(1).FaceColor=blue;bb(2).FaceColor=orange;yline(0,'-','Color',[.2 .2 .2]);yline(10,':','Color',gray);yline(-10,':','Color',gray);xticks(1:4);xticklabels(piles);ylabel('Signed prediction error (%)');panel(k);if k==1,legend({'Four reactions','{\it p} only'},'Location','northeast');end;end
saveall(f,'Fig23_field_metric_errors',2,'Field errors without parameter selection. (a) Initial secant stiffness at ground displacement D/1000. (b) Load at D/10. Signed error is prediction/observation minus one. The four-reaction and p-only models share geometry, input profiles, beam theory and loading conditions; p-only coefficients are not recalibrated. Dotted limits show the 10% research screening criterion. The complete model improves some metrics but is not uniformly more accurate; in particular, its larger Cowden load at D/10 increases overprediction. Four specimen results remain grouped into two site pairs.','field_metrics.csv','Shows both improvement and deterioration when adding calibrated reaction channels.');

f=canvas(2,2,18,13);
for site=1:2
 for metric=1:2
  nexttile((site-1)*2+metric);hold on;
  if site==1,ps=["DL1","DL2"];sc=["API_phi32","API_phi36","API_phi40"];xx=[32 36 40];xl='Assumed peak friction angle, {\it\phi} (deg)';else,ps=["CL1","CL2"];sc=["Matlock_e004_J050","Matlock_e005_J050","Matlock_e007_J050"];xx=[.004 .005 .007];xl='Assumed strain parameter, {\it\epsilon}_{50}';end
  for j=1:2,v=zeros(1,3);for k=1:3,q=E(string(E.pile_id)==ps(j)&string(E.scenario_id)==sc(k),:);if metric==1,v(k)=q.Kinit_relative_error;else,v(k)=q.H01D_relative_error;end;end;plot(xx,100*v,'-o','Color',cols(j,:),'MarkerSize',4);end
  yline(0,':k');xlabel(xl);ylabel('Signed prediction error (%)');if site==1&&metric==2,legend(ps,'Location','southeast');else,legend(ps,'Location','northeast');end;xlim([min(xx)-.03*range(xx),max(xx)+.03*range(xx)]);panel((site-1)*2+metric);
 end
end
saveall(f,'Fig24_baseline_parameter_dependence',4,'Dependence of conventional baseline predictions on unmeasured constitutive inputs. (a,b) Dunkirk initial-stiffness and D/10-load errors for assumed peak friction angles of 32, 36 and 40 degrees, using the documented API stiffness-gradient relation. The measured critical-state angle is not treated as a measured peak angle. (c,d) Corresponding Cowden errors for epsilon50=0.004, 0.005 and 0.007 at J=0.5. These strain values are generic overconsolidated-clay scenarios; the Matlock soft-clay framework remains an extrapolative baseline for stiff glacial till. Both piles are shown, and no scenario is chosen because it produces the best response match.','field_metrics.csv; parameter_scenarios.csv','Exposes missing-input confounding and the possibility of matching capacity while missing stiffness.');

f=canvas(2,2,18,13);
for j=1:4
 nexttile;hold on;plot([-10 10 10 -10 -10],[-10 -10 10 10 -10],':','Color',gray);plot(0,0,'k+','MarkerSize',8,'LineWidth',1.2);
 q=E(string(E.pile_id)==piles(j)&string(E.scenario_id)=="PISA_full",:);plot(100*q.Kinit_relative_error,100*q.H01D_relative_error,'ko','MarkerFaceColor','k','MarkerSize',5);
 for family=1:2
  if family==1,ss=["PISA_G0_080","PISA_full","PISA_G0_120"];else,ss=["PISA_Su_080","PISA_full","PISA_Su_120"];if j<=2,continue;end;end
  vals=[];for k=1:3,q=E(string(E.pile_id)==piles(j)&string(E.scenario_id)==ss(k),:);vals=[vals;100*q.Kinit_relative_error,100*q.H01D_relative_error];end
  plot(vals(:,1),vals(:,2),'-s','Color',cols(family,:),'MarkerSize',4);
 end
 xlabel('Initial-stiffness error (%)');ylabel('{\it D}/10 load error (%)');panel(j);xlim([-40 20]);ylim([-30 45]);if j==3,legend({'10% screening window','Published metrics','Nominal','Modulus scenarios','Strength scenarios'},'Location','northwest','FontSize',7.5);end
end
saveall(f,'Fig25_field_scenario_tradeoffs',4,'Joint stiffness and load discrepancies across predeclared PISA input scenarios. Panels (a–d) show DL1, DL2, CL1 and CL2. Black circles are nominal predictions; blue paths vary G0 by factors 0.8, 1.0 and 1.2, while orange paths vary Cowden su by the same factors. The origin corresponds to the published two-metric match and the square is a 10% research screening window. Scenario paths are not confidence regions. Failure to reach the origin within these specified perturbations is reported without changing the nominal inputs or fitting the test curves.','field_metrics.csv','Tests whether simple stiffness or strength adjustments explain both discrepancies simultaneously.');

f=canvas(2,2,18,13);
for site=1:2
 if site==1,ps=["DL1","DL2"];else,ps=["CL1","CL2"];end
 for metric=1:2
  nexttile((metric-1)*2+site);hold on;
  for j=1:2
   q=P(string(P.pile_id)==ps(j)&string(P.scenario_id)=="PISA_full",:);if metric==1,xx=q.p_N_per_m/1e6;else,xx=abs(q.M_theta_Nm)/1e6;end;plot(xx,q.z_m,'Color',cols(j,:),'LineStyle',char(repmat('-',1,j)));
  end
  set(gca,'YDir','reverse');ylabel('Depth, {\it z} (m)');if metric==1,xlabel('Distributed force, {\it p} (MN/m)');else,xlabel('Bending moment, {\it M} (MN m)');end;panel((metric-1)*2+site);if metric==1,legend(ps,'Location','west');end
 end
end
saveall(f,'Fig26_field_load_transfer',4,'Predicted internal load transfer at common ground displacement D/10. Columns show the Dunkirk and Cowden pairs. (a,b) Signed distributed lateral reaction. (c,d) Bending-moment magnitude. Each pile retains its actual embedment, wall thickness and loading height; the resulting head forces therefore differ. These profiles are model predictions and have not been matched to measured internal soil reactions or structural strain profiles.','ultimate_profiles.csv','Connects differences in head response to predicted internal force distributions without presenting them as measurements.');

writetable(cell2table(reg,'VariableNames',{'id','panels','caption','data_sources','evidence_value','program'}),fullfile(out,'field_figure_captions_evidence.csv'));
 function f=canvas(nr,nc,w,h),f=figure('Visible','off','Units','centimeters','Position',[1 1 w h],'Color','w','Renderer','painters');tiledlayout(nr,nc,'TileSpacing','compact','Padding','compact');end
 function panel(k),title('');text(.015,.985,['(' char('a'+k-1) ')'],'Units','normalized','VerticalAlignment','top','FontSize',10,'FontWeight','bold','Tag','PanelLabel');end
 function saveall(f,id,np,caption,data,purpose)
  for ax=reshape(findall(f,'Type','axes'),1,[]),ax.FontSize=9;ax.XLabel.FontSize=10;ax.YLabel.FontSize=10;end
  for lg=reshape(findall(f,'Type','legend'),1,[]),lg.FontSize=8.5;end
  for tx=reshape(findall(f,'Tag','PanelLabel'),1,[]),tx.FontSize=10.5;end
  drawnow;for ax=reshape(findall(f,'Type','axes'),1,[]),ax.Box='off';end
  savefig(f,fullfile(out,[id '.fig']));exportgraphics(f,fullfile(out,[id '.pdf']),'ContentType','vector');print(f,fullfile(out,[id '.svg']),'-dsvg');print(f,fullfile(out,[id '.emf']),'-dmeta');exportgraphics(f,fullfile(out,[id '_QA.png']),'Resolution',180);reg(end+1,:)={id,np,caption,data,purpose,'code/plot_field_figures.m'};close(f);
 end
end
