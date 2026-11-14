function out=plot_evaluation_results(auditDir)

here=fileparts(mfilename('fullpath'));
if nargin<1,auditDir=fullfile(here,'reference_results');end
out=fullfile(auditDir,'figures');if ~isfolder(out),mkdir(out);end
T=readtable(fullfile(auditDir,'null_test_prediction_and_interval.csv'));
A=readtable(fullfile(auditDir,'independent_query_and_total_cost_audit.csv'),'TextType','string');
S=readtable(fullfile(auditDir,'cost_and_quality_summary.csv'),'TextType','string');
E=readtable(fullfile(auditDir,'paired_adaptive_minus_periodic_cost.csv'));
record=jsondecode(fileread(fullfile(auditDir,'audit_record.json')));n=record.null_control;
set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman');

rc=[0 .447 .698;.9 .48 0;0 .62 .45];rm={'o','s','^'};
f=canvas();t=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
a=nexttile(t);hold(a,'on');h=gobjects(3,1);
for j=1:3,r=T.constraint==j;h(j)=plot(a,T.direct_prediction(r),T.source_prediction(r),rm{j},'Color',rc(j,:),'MarkerSize',3,'LineWidth',.65);end
v=[min(T.direct_prediction) max(T.direct_prediction)];plot(a,v,v,'k--','LineWidth',.7);xlabel(a,'Direct reference prediction');ylabel(a,'Source-corrected prediction');
legend(a,h,{'Rotation','Displacement','Stress'},'Location','north','Box','off','FontSize',7.5);panel(a,'a');
b=nexttile(t);hold(b,'on');
for j=1:3,r=T.constraint==j;semilogy(b,T.row(r),max(abs(T.source_minus_direct(r)),1e-16),rm{j},'Color',rc(j,:),'MarkerSize',2.6,'LineWidth',.6);end
set(b,'YScale','log');xlabel(b,'Held-out row index');ylabel(b,'Absolute prediction difference');yline(b,1e-10,'k:','LineWidth',.7);ylim(b,[1e-16 1e-9]);panel(b,'b');
c=nexttile(t);values=[n.max_operator_additivity_absolute_error,n.max_reference_operator_implementation_error,n.max_saved_source_minus_direct_absolute_error];
plot(c,1:3,max(values,1e-16),'ko','MarkerFaceColor',[.45 .25 .65],'MarkerSize',4);set(c,'YScale','log','XTick',1:3,'XTickLabel',{'Additivity','Direct vs. {\bf\it S}{\bf\it Y}','Saved means'});xlim(c,[.5 3.5]);ylim(c,[1e-16 1e-9]);yline(c,1e-10,'k:','LineWidth',.7);ylabel(c,'Maximum absolute residual');panel(c,'c');
d=nexttile(t);ratios=zeros(3,1);for j=1:3,ratios(j)=mean(T.width_ratio(T.constraint==j));end
hold(d,'on');for j=1:3,bar(d,j,ratios(j),.8,'FaceColor',rc(j,:),'EdgeColor',rc(j,:));end
yline(d,1,'k:','LineWidth',.7);xticks(d,1:3);xticklabels(d,{'Rotation','Displacement','Stress'});ylabel(d,'Source / direct interval half-width');ylim(d,[0 max(1.25,max(ratios)*1.12)]);panel(d,'d');
formatExport(f,out,'E01_additive_smoother_null');

mc=[.3 .3 .35;0 .447 .698;.9 .48 0;0 .62 .45;.45 .25 .65];
marks={'x','s','o','d','^'};styles={'--','--','-','--','-'};
names=["Frozen","LowOnly","ReferenceOnly","Periodic","DecisionAdaptive"];
short={'Frozen','Low','Reference','Periodic','Adaptive'};
f=canvas();t=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
a=nexttile(t);s=S(S.query_budget==16,:);values=[s.mean_initial_seconds,s.mean_online_seconds,s.mean_uniform_audit_seconds,s.mean_fit_seconds];
bh=bar(a,values,'stacked','EdgeColor','none');cc=[.78 .78 .78;.4 .4 .4;.92 .78 .48;.45 .55 .67];for k=1:4,bh(k).FaceColor=cc(k,:);end
xticks(a,1:5);xticklabels(a,short);xtickangle(a,25);ylabel(a,'Solver estimate + measured fit time (s)');ylim(a,[0 max(sum(values,2))*1.08]);
legend(a,bh,{'Initial','Online','Final audit','Fitting'},'Location','southoutside','NumColumns',2,'Box','off','FontSize',7.5);panel(a,'a');
b=nexttile(t);hold(b,'on');h=gobjects(5,1);
for k=1:5,r=string(A.method)==names(k)&A.query_budget==16;h(k)=plot(b,A.posthoc_solver_plus_fit_seconds(r),100*A.conditional_regret(r),marks{k},'Color',mc(k,:),'MarkerSize',4,'LineWidth',.8);end
xlabel(b,'Post hoc solver estimate + fits (s)');ylabel(b,'Feasible mass gap (%)');ylim(b,[-2 55]);
legend(b,h,short,'Location','north','Box','off','FontSize',7.5);panel(b,'b');
c=nexttile(t);hold(c,'on');
for k=1:5,s=S(string(S.method)==names(k),:);[~,ii]=sort(s.query_budget);s=s(ii,:);plot(c,s.mean_posthoc_total_seconds,100*s.mean_conditional_regret,'LineStyle',styles{k},'Marker',marks{k},'Color',mc(k,:),'MarkerSize',4,'LineWidth',.9);end
xlabel(c,'Mean post hoc cost (s)');ylabel(c,'Mean feasible mass gap (%)');ylim(c,[-2 55]);panel(c,'c');
d=nexttile(t);hold(d,'on');bud=unique(E.query_budget);
for i=1:numel(bud),v=E.adaptive_minus_periodic_posthoc_total_seconds(E.query_budget==bud(i));plot(d,i+linspace(-.12,.12,numel(v))',v,'^','Color',mc(5,:),'MarkerSize',4,'LineWidth',.8);plot(d,i,mean(v),'kd','MarkerFaceColor','k','MarkerSize',5);end
yline(d,0,'k:','LineWidth',.7);xlim(d,[.5 numel(bud)+.5]);ylim(d,[-.005 max(E.adaptive_minus_periodic_posthoc_total_seconds)*1.12]);xticks(d,1:numel(bud));xticklabels(d,string(bud));xlabel(d,'Additional query budget (R equivalents)');ylabel(d,'Adaptive minus periodic cost (s)');panel(d,'d');
formatExport(f,out,'E02_query_and_total_cost');

captions={
'E01_additive_smoother_null','H2','All','Additive source correction as a linear-smoother null control. All three models use the same 80 training geometries, Gaussian-kernel length scale 2.4 and ridge 0.001. (a) Source-corrected versus direct-reference predictions of rotation, displacement and stress utilization for 120 independent test geometries; the dashed line denotes equality. (b) Absolute prediction differences by test-row index. Values below 10⁻¹⁶ are displayed at that floor; dotted lines in (b,c) mark the 10⁻¹⁰ numerical tolerance. (c) Maximum residuals of the additive operator identity, explicit operator predictions versus saved direct predictions, and saved source versus direct predictions, evaluated over all 260 geometry coordinates and three responses. (d) Mean source/direct calibrated half-width ratios by response on the test set; the dotted line denotes equal width. This post-evaluation diagnostic uses frozen predictions without refitting. Equal point predictions arise from the shared linear smoother, whereas joint calibrated intervals can differ because their response-wise scale ratios differ.','Post-evaluation mathematical diagnostic','null_test_prediction_and_interval.csv; null_model_parameters.csv','SI, or main if the source-decomposition claim is discussed';
'E01_additive_smoother_null','H2','a','Source-corrected versus direct-reference predictions for rotation, mudline displacement and stress utilization. The dashed line is equality. Colours and marker shapes identify responses.','Shows pointwise equality rather than coincident aggregate RMSE','null_test_prediction_and_interval.csv','SI';
'E01_additive_smoother_null','H2','b','Absolute pointwise prediction differences for the 120 held-out geometry rows. Values below 1e-16 are displayed at 1e-16 only to permit a logarithmic axis. The dotted line is the 1e-10 numerical-check tolerance.','Checks response-wise equality across all held-out cases','null_test_prediction_and_interval.csv','SI';
'E01_additive_smoother_null','H2','c','Maximum absolute residuals for S Y_L + S(Y_R-Y_L) - S Y_R, for explicit S Y_R versus saved direct predictions, and for saved source versus direct means. The maximum is over all 260 available geometry coordinates and three responses. The dotted line is the 1e-10 tolerance.','Separates operator identity from floating-point implementation consistency','audit_record.json','SI';
'E01_additive_smoother_null','H2','d','Mean source/direct calibrated half-width ratio by response on the held-out geometries. The dotted line is equal width. A common residual field does not require equal joint intervals because the response-wise scale ratios differ.','Explains why equal means do not imply equal joint decisions','null_test_prediction_and_interval.csv','SI';
'E02_query_and_total_cost','H3','All','Physical-query accounting and computational cost on the 128-candidate numerical domain. Solver estimates use independent development medians of 0.00354415 s per low query and 0.02267415 s per reference query; fitting times were measured. (a) Mean initialization, online-query, final-audit and fitting costs across eight shared initializations at additional query budget 16. Initialization comprises 12 low and four reference queries. The uniform final audit adds one logical reference fee, including previously queried selections. (b) Individual cost–mass-gap outcomes at the same query cap; mass gaps are relative to the candidate-domain reference optimum. (c) Mean cost–quality outcomes at caps 8, 16 and 24; lines connect budget conditions. (d) Paired adaptive-minus-periodic costs: triangles show eight initializations per cap, black diamonds their means, and the dotted line equality. Horizontal offsets separate points. All 120 final recommendations are feasible under the numerical reference. These post hoc solver-plus-fit totals exclude offline candidate-table generation and are not measured end-to-end execution times. A common solver-query budget does not imply equal total computation.','Post-evaluation cost decomposition; not a new efficiency-optimized algorithm','independent_query_and_total_cost_audit.csv; cost_and_quality_summary.csv; paired_adaptive_minus_periodic_cost.csv','SI';
'E02_query_and_total_cost','H3','a','Mean cost components at additional query budget 16 across eight shared initializations: initial queries, online queries, uniform final audit and measured fitting. Component colours encode cost categories and are independent of the policy colour map.','Quantifies fitting outside the solver-query budget','cost_and_quality_summary.csv','SI';
'E02_query_and_total_cost','H3','b','Cost and feasible mass gap for all eight initializations of every policy at additional query budget 16. The mass gap is relative to the optimum within the 128-candidate numerical domain. Policy colours are the same in panels b and c.','Displays individual outcomes at matched solver-query caps','independent_query_and_total_cost_audit.csv','SI';
'E02_query_and_total_cost','H3','c','Mean estimated cost and mean feasible mass gap for additional query caps 8, 16 and 24. Lines join the three budget conditions of the same policy; they are not fitted continuous performance curves. The Frozen outcomes overlap because no online queries are performed.','Shows cost-quality trade-offs without claiming matched total computation','cost_and_quality_summary.csv','SI';
'E02_query_and_total_cost','H3','d','Paired adaptive-minus-periodic post hoc costs at each budget. Open triangles show the eight shared initializations, with horizontal offsets for readability; black diamonds show their means. The dotted line is equal cost. All 24 paired differences are positive in this run.','Supplies a counterexample to interpreting equal query budgets as equal or superior total computational cost','paired_adaptive_minus_periodic_cost.csv','SI'};
C=cell2table(captions,'VariableNames',{'figure','hypothesis','panel','caption','evidence_value','data_sources','placement'});writetable(C,fullfile(out,'figure_captions_and_evidence.csv'));
style=struct('font','Times New Roman','ordinary_pt',9,'legend_pt',7.5,'panel_pt',10,'panel_fontweight','bold','panel_units','normalized','panel_xy',[.015 .985],...
 'titles',false,'legends_box','off','canvas_cm',[18.4 14.6],'response_colors',rc,'policy_colors',mc,'raster_images_are_QA_only',true);
fid=fopen(fullfile(out,'figure_style.json'),'w','n','UTF-8');cl=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(style,'PrettyPrint',true),'char');
disp(out);
end
function f=canvas(),f=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 18.4 14.6]);end
function panel(a,letter)
text(a,.015,.985,['(' letter ')'],'Units','normalized','FontName','Times New Roman','FontSize',10,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','top','Clipping','off','Tag','PanelLabel');
end
function formatExport(f,out,name)
a=findall(f,'Type','axes');set(a,'FontName','Times New Roman','FontSize',9,'FontAngle','normal','TickLabelInterpreter','tex','LineWidth',.65,'TickDir','out','Box','off');tx=findall(f,'Type','text');set(tx,'FontName','Times New Roman','FontAngle','normal');set(findall(f,'Type','legend'),'Interpreter','tex','FontAngle','normal');set(f,'Renderer','painters');
labels=findall(f,'Tag','PanelLabel');set(labels,'FontSize',10,'FontWeight','bold');drawnow;p=fullfile(out,name);savefig(f,[p '.fig']);print(f,[p '.svg'],'-dsvg');exportgraphics(f,[p '.pdf'],'ContentType','vector');print(f,[p '.emf'],'-dmeta');exportgraphics(f,[p '_QA.png'],'Resolution',180);close(f);
end
