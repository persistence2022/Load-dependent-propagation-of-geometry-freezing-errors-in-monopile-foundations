function plot_auxiliary_results(runDir)

root=fileparts(mfilename('fullpath'));addpath(root);
out=fullfile(runDir,'figures');if ~isfolder(out),mkdir(out);end
S=load(fullfile(runDir,'data_formal.mat'));d=S.data;c=analysis_config();
set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman','defaultAxesFontSize',9,...
 'defaultTextFontSize',9,'defaultAxesLineWidth',.7,'defaultLineLineWidth',1.3,'defaultFigureColor','w','defaultAxesTickDir','out');
colors=[0 .447 .698;.9 .48 .0;0 .62 .45;.45 .25 .65];labels={'Frozen secant','Updated secant','Frozen nonlinear','Updated nonlinear'};modelStyles={'--','--','-','-'};
captions={};specs={};

f=canvas(18,12);tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;hold on;role=["dev","cal","test"];
for k=1:3,ii=d.split==role(k);scatter(d.X(ii,1),d.X(ii,2)*1e3,14,colors(k,:),'filled','MarkerFaceAlpha',.65);end
xlabel('Diameter (m)');ylabel('Wall thickness (mm)');ylim([60 150]);yticks(60:20:140);legend({'Development','Calibration','Independent test'},'Location','best','FontSize',7);ttl('a  Geometry partitions');
nexttile;scatter(d.X(:,3),d.mass/1e3,18,max(d.reference,[],2),'filled');xlabel('Embedded length (m)');ylabel('Embedded steel mass (t)');cb=colorbar;cb.Label.String='Reference maximum utilization';ttl('b  Engineering range');
nexttile;hold on;for j=1:3,[x,y]=ecdf0(d.reference(:,j));plot(x,y,'Color',colors(j,:));end;xline(1,':k');xlabel('Reference utilization');ylabel('Empirical cumulative fraction');legend({'Rotation','Displacement','Stress'},'Location','southeast','FontSize',7);ttl('c  Constraint coverage');
nexttile;hold on;costbars=median(d.elapsed,1)*1e3;for a=1:4,bar(a,costbars(a),.8,'FaceColor',colors(a,:));end;xticks(1:4);xticklabels({'FS','US','FN','UN'});ylabel('Median paired-load solve time (ms)');ttl('d  Measured model cost');grid on;
savefigs(f,'Fig02_design_and_cost');cap('Fig02_design_and_cost','H1;H2;H3','Independent geometry partitions and model coverage. (a) Development, calibration and held-out test geometries. (b) Embedded steel mass and reference constraint utilization. (c) Cumulative utilization distributions for the three declared constraints. (d) Actual median solver times for both load cases per geometry; FS/US denote frozen/updated half-capacity secant models and FN/UN the corresponding nonlinear models. These are controlled numerical scenarios.','data_formal.csv','Each panel identifies a different design, constraint or cost dimension; paired load cases are not independent observations.');
E=load(fullfile(runDir,'evaluation','evaluation_result.mat'));e=E.result;nm={'Frozen low','Global correction','Direct reference','Source correction'};
f=canvas(18,6.6);tiledlayout(1,3,'TileSpacing','compact','Padding','compact');C=e.curves;P=e.paired_effects;
nexttile;hold on;keys={'low_proxy','global_p95','direct_reference','source_corrected'};for a=1:4,ii=string(C.method)==keys{a};plot(C.approved_fraction(ii),C.false_approval_given_approved(ii),'-o','MarkerSize',3,'Color',colors(a,:));end;xlabel('Approved fraction');ylabel('False approvals / approved');legend(nm,'FontSize',7,'Location','northwest');ttl('a  Matched acceptance');
nexttile;hold on;for a=1:4,keys={'low_proxy','global_p95','direct_reference','source_corrected'};ii=string(C.method)==keys{a};plot(C.approved_fraction(ii),C.mean_positive_violation_approved(ii),'-o','MarkerSize',3,'Color',colors(a,:));end;xlabel('Approved fraction');ylabel('Mean positive violation');ttl('b  Engineering severity');
nexttile;ii=string(P.comparison)=="source_corrected_minus_global_p95";errorbar(P.target_fraction(ii),P.false_approval_difference(ii),P.false_approval_difference(ii)-P.difference_ci_low(ii),P.difference_ci_high(ii)-P.false_approval_difference(ii),'-o','Color',colors(4,:),'MarkerSize',3);yline(0,':k');xline(.5,':k');ylim([min(P.difference_ci_low(ii))-.005 .005]);xlabel('Approved fraction');ylabel('Source minus global risk');ttl('c  Paired group uncertainty');
savefigs(f,'Fig09_matched_decision_quality');cap('Fig09_matched_decision_quality','H2','Decision quality at matched approval fractions. (a) False approvals among the accepted designs. (b) Mean positive maximum-constraint violation among accepted designs. (c) Paired source-corrected minus global-correction differences with 95% percentile intervals from geometry-group resampling. The 50% approval comparison was specified before formal evaluation; the full curve is secondary. Forced acceptance for ranking comparison is distinct from interval-based approval.','evaluation/matched_approval_curves.csv; evaluation/paired_effects.csv','Equal acceptance avoids an apparent benefit caused only by withholding more decisions.');

P=load(fullfile(runDir,'policy','policy_result.mat'));r=P.result;T=r.table;SS=r.summary;methods=string(c.policy.methods);mc=[.3 .3 .35;0 .447 .698;.9 .48 0;0 .62 .45;.45 .25 .65];
f=canvas(18,12);tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;hold on;for a=1:5,frac=zeros(1,3);for b=1:3,ii=string(T.method)==methods(a)&T.additional_budget_R_equivalent==c.policy.budget_reference_equivalents(b);frac(b)=mean(T.feasible_mass_regret_fraction(ii)<1e-8);end;plot(c.policy.budget_reference_equivalents,frac,'-o','Color',mc(a,:),'MarkerSize',4);end;ylim([0 1.05]);xlabel('Additional solver-query budget (R calls)');ylabel('Fraction reaching candidate optimum');legend({'Frozen','Low queries','Reference queries','Periodic','Adaptive'},'FontSize',7,'Location','best');ttl('a  Candidate-optimum attainment');
nexttile;hold on;for a=1:5,ii=string(SS.method)==methods(a);plot(SS.budget(ii),100*SS.mean_conditional_regret(ii),'-o','Color',mc(a,:),'MarkerSize',4);end;xlabel('Additional solver-query budget (R calls)');ylabel('Feasible mass gap (%)');ttl('b  Conditional material outcome');
nexttile;hold on;for a=1:5,ii=string(T.method)==methods(a)&T.additional_budget_R_equivalent==16;scatter(repmat(a,sum(ii),1)+linspace(-.12,.12,sum(ii))',100*T.feasible_mass_regret_fraction(ii),24,mc(a,:),'filled');end;xticks(1:5);xticklabels({'Frozen','Low','Ref.','Periodic','Adaptive'});ylabel('Feasible mass gap (%)');ttl('c  All initializations, middle budget');
nexttile;hold on;for a=1:5,ii=string(SS.method)==methods(a);plot(SS.mean_total_query_cost(ii),SS.mean_fit_seconds(ii),'-o','Color',mc(a,:),'MarkerSize',4);end;xlabel('Total accounted query cost (R calls)');ylabel('Model-fitting time (s)');ttl('d  Fitting cost outside query budget');
savefigs(f,'Fig10_budgeted_design_outcomes');cap('Fig10_budgeted_design_outcomes','H3','Finite-candidate design outcomes for five query policies and eight shared initializations. (a) Candidate-optimum attainment using a relative mass tolerance of 1e-8. (b) Mean feasible mass gap. (c) Individual outcomes at the middle budget; horizontal offsets separate repeated runs without changing their values. (d) Measured model-fitting time. All 120 final recommendations satisfy the numerical reference constraints. Budgets constrain solver-query equivalents using measured development costs, not total computation including refitting. The common final reference audit is accounted separately. Repeated initializations do not represent independent sites.','policy/policy_results.csv; policy/policy_summary.csv','Material outcomes and separate query and fitting costs.');

writetable(cell2table(captions,'VariableNames',{'figure','hypotheses','caption','data_sources','evidence_value'}),fullfile(out,'figure_captions_and_evidence.csv'));
copyfile(fullfile(runDir,'data_formal.csv'),fullfile(out,'Fig02_Fig03_source_data.csv'));
stage.save_json(fullfile(out,'figure_style.json'),struct('font','Times New Roman','axes_font_pt',9,'canvas_cm',[18 12],'colors',colors,'native_formats',{{'FIG','SVG','PDF','EMF'}},'raster_previews_only',true));
 function f=canvas(w,h),f=figure('Visible','off','Units','centimeters','Position',[1 1 w h],'Renderer','painters');end
 function ttl(s),title('');text(.015,.985,['(' s(1) ')'],'Units','normalized','VerticalAlignment','top','FontWeight','bold','FontSize',10,'Interpreter','none');box off;end
 function savefigs(f,name)
  if startsWith(name,'Fig08')||startsWith(name,'Fig09')
   lsh={'--',':','-','--'};mkh={'s','d','o','^'};ll=findall(f,'Type','line');
   for obj=reshape(ll,1,[]),[distance,idx]=min(sum((colors-obj.Color).^2,2));if distance<1e-10,obj.LineStyle=lsh{idx};obj.Marker=mkh{idx};obj.MarkerIndices=unique(round(linspace(1,numel(obj.XData),min(8,numel(obj.XData)))));obj.MarkerSize=3.5;end,end
  end
  if startsWith(name,'Fig10')||startsWith(name,'Fig11')
   ls0={'--','--','-','--','-'};mk0={'x','s','o','d','^'};lines=findall(f,'Type','line');
   for obj=reshape(lines,1,[]),[distance,idx]=min(sum((mc-obj.Color).^2,2));if distance<1e-10,obj.LineStyle=ls0{idx};obj.Marker=mk0{idx};obj.MarkerIndices=unique(round(linspace(1,numel(obj.XData),min(7,numel(obj.XData)))));obj.MarkerSize=4;end,end
  end
  ls=findall(f,'Type','legend');for lg=reshape(ls,1,[]),lg.Box='off';lg.FontSize=7.5;end
  if strcmp(name,'Fig02_design_and_cost')
   for lg=reshape(ls,1,[]),if any(contains(string(lg.String),'Development')),lg.Orientation='horizontal';lg.NumColumns=3;lg.Layout.Tile='north';end,end
  elseif strcmp(name,'Fig10_budgeted_design_outcomes')
   ls.Orientation='horizontal';ls.NumColumns=5;ls.Layout.Tile='north';
   aa=findall(f,'Type','axes');for ax=reshape(aa,1,[]),if any(strcmp(string(ax.XTickLabel),'Periodic')),ax.XTickLabelRotation=20;end,end
  end
  drawnow;stage.vector_colorbars(f);drawnow;
  fullCanvas=ismember(name,{'Fig02_design_and_cost','Fig10_budgeted_design_outcomes'});
  if fullCanvas,f.PaperUnits='centimeters';f.Units='centimeters';fp=f.Position;f.PaperPosition=[0 0 fp(3:4)];f.PaperSize=fp(3:4);f.PaperPositionMode='manual';annotation(f,'rectangle',[0.002 0.002 0.996 0.996],'Color',[1 1 1],'LineWidth',0.1,'Tag','PublicationCanvasBoundary');drawnow;end
  savefig(f,fullfile(out,[name '.fig']));exportgraphics(f,fullfile(out,[name '.pdf']),'ContentType','vector','BackgroundColor','white');print(f,fullfile(out,[name '.svg']),'-dsvg');print(f,fullfile(out,[name '.emf']),'-dmeta');exportgraphics(f,fullfile(out,[name '_QA.png']),'Resolution',180);close(f);
 end
 function cap(id,h,text,src,value),captions(end+1,:)={id,h,text,src,value};end
end
function [x,y]=ecdf0(x),x=sort(x(isfinite(x)));y=(1:numel(x))'/max(numel(x),1);end
function m=divmap(),a=[0 .447 .698];b=[.95 .95 .93];c=[.9 .42 .02];m=[a+(b-a).*linspace(0,1,128)';b+(c-b).*linspace(0,1,128)'];end
function drawcells(x,y,Z,map)

dx=mean(diff(x));dy=mean(diff(y));hold on;
for i=1:numel(x),for j=1:numel(y),rectangle('Position',[x(i)-dx/2,y(j)-dy/2,dx,dy],'FaceColor',map(Z(j,i),:),'EdgeColor','none');end,end
xlim([min(x)-dx/2 max(x)+dx/2]);ylim([min(y)-dy/2 max(y)+dy/2]);colormap(gca,map);clim([.5 4.5]);box on;
end
