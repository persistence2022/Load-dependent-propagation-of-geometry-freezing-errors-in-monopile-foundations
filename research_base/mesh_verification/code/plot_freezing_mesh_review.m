function plot_freezing_mesh_review

review=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(review,'figures');if ~isfolder(out),mkdir(out);end
diary(fullfile(out,'Fig36_plot_matlab.log'));cleanup=onCleanup(@()diary('off'));
fprintf('Plot-only execution, MATLAB %s, %s\n',version,computer);
R=readtable(fullfile(review,'summary','targeted_mesh_responses_combined.csv'),'TextType','string','Delimiter',',');
B=readtable(fullfile(review,'summary','targeted_freezing_bias_combined.csv'),'TextType','string','Delimiter',',');
F=readtable(fullfile(review,'reference_results','strict','failed_cases.csv'),'TextType','string','Delimiter',',');
fprintf('Input rows, responses %d, pairs %d, retained failures %d\n',height(R),height(B),height(F));
assert(height(R)==8&&height(B)==4&&height(F)==2);
cases=["Clay_D1","Clay_D2","Sand_D2t_Dr55","Sand_D2t_Dr85"];
colors=[0 114 178;213 94 0;0 158 115;204 121 167]/255;
P=table;Hheight=[37.5 87.5 87.5 87.5];
for j=1:4
    u=R(R.case_id==cases(j)&R.treatment=="Updated",:);
    f=R(R.case_id==cases(j)&R.treatment=="Frozen",:);
    b=B(B.case_id==cases(j),:);assert(height(u)==1&&height(f)==1&&height(b)==1);
    b40=100*(f.v40_m/u.v40_m-1);b80=100*(f.v80_m/u.v80_m-1);
    du=100*(u.v80_m/u.v40_m-1);df=100*(f.v80_m/f.v40_m-1);
    assert(abs(b40-b.bias40_percent)<1e-10&&abs(b80-b.bias80_percent)<1e-10);
    assert(abs((b80-b40)-b.bias_change_percentage_points)<1e-10);
    assert(abs(du-100*u.relative_displacement_change)<1e-10&&abs(df-100*f.relative_displacement_change)<1e-10);
    failed=sum(F.case_id==cases(j));
    row=table(cases(j),u.D_m,u.H_N,Hheight(j)/u.D_m,...
        u.v40_m,f.v40_m,u.v80_m,f.v80_m,b40,b80,b80-b40,du,df,...
        u.actual_geometry_domain,u.diagnostic_tolerance,u.max_relative_residual,f.max_relative_residual,...
        failed,double(failed>0),u.v80_m<=.1*u.D_m,f.v80_m<=.1*f.D_m,...
        'VariableNames',{'case_id','D_m','H_N','h_over_D','updated_v40_m','frozen_v40_m',...
        'updated_v80_m','frozen_v80_m','bias40_percent','bias80_percent','bias_change_pp',...
        'updated_displacement_change_percent','frozen_displacement_change_percent',...
        'actual_geometry_domain','accepted_80_element_residual_tolerance','updated_max_relative_residual',...
        'frozen_max_relative_residual','retained_strict_failed_branches','diagnostic_retry',...
        'updated_80_within_D_over_10','frozen_80_within_D_over_10'});
    P=[P;row];
end
assert(all(abs(P.bias_change_pp)<=.5));
assert(all(abs([P.updated_displacement_change_percent;P.frozen_displacement_change_percent])<=1));
writetable(P,fullfile(out,'Fig36_mesh_effect_on_freezing_plot_data.csv'));
writetable(R,fullfile(out,'Fig36_complete_response_data.csv'));
writetable(F,fullfile(out,'Fig36_retained_strict_failures.csv'));

set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman',...
    'defaultAxesFontSize',9,'defaultTextFontSize',9,'defaultAxesLineWidth',.7,...
    'defaultLineLineWidth',1.25,'defaultAxesTickDir','out','defaultAxesBox','off',...
    'defaultLegendBox','off','defaultLegendFontSize',8.5);
fig=figure('Visible','off','Color','w','Units','centimeters',...
    'Position',[1 1 18 10],'Renderer','painters');
tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
tl.OuterPosition=[0 .20 1 .80];
ticks={'Clay\newlineD1*','Clay\newlineD2','Sand\newline55%','Sand\newline85%'};

ax=nexttile;hold on;
for j=1:4
    plot([j-.10,j+.10],[P.bias40_percent(j),P.bias80_percent(j)],'-','Color',colors(j,:));
    plot(j-.10,P.bias40_percent(j),'o','Color',colors(j,:),'MarkerFaceColor',colors(j,:),'MarkerSize',4.5);
    plot(j+.10,P.bias80_percent(j),'s','Color',colors(j,:),'MarkerFaceColor','w','MarkerSize',4.5);
end
a=plot(nan,nan,'ko','MarkerFaceColor','k','MarkerSize',4.5);
b=plot(nan,nan,'ks','MarkerFaceColor','w','MarkerSize',4.5);
legend([a b],{'40 elements','80 elements'},'Location','southwest');
ylabel({'Frozen-diameter','displacement bias (%)'});ylim([18 66]);
decorate(ax,ticks,'a');

ax=nexttile;hold on;yline(0,':','Color',[.5 .5 .5]);
for j=1:4
    plot([j-.10,j+.10],[P.updated_displacement_change_percent(j),P.frozen_displacement_change_percent(j)],'-','Color',colors(j,:));
    plot(j-.10,P.updated_displacement_change_percent(j),'o','Color',colors(j,:),'MarkerFaceColor',colors(j,:),'MarkerSize',4.5);
    plot(j+.10,P.frozen_displacement_change_percent(j),'s','Color',colors(j,:),'MarkerFaceColor','w','MarkerSize',4.5);
end
a=plot(nan,nan,'ko','MarkerFaceColor','k','MarkerSize',4.5);
b=plot(nan,nan,'ks','MarkerFaceColor','w','MarkerSize',4.5);
legend([a b],{'Updated','Frozen'},'Location','northeast');
ylabel({'Displacement change','from 40 to 80 elements (%)'});ylim([-.014 .046]);
ax.YAxis.Exponent=0;ytickformat('%.2f');
decorate(ax,ticks,'b');

ax=nexttile;hold on;yline(0,':','Color',[.5 .5 .5]);
for j=1:4
    plot([j j],[0 P.bias_change_pp(j)],'-','Color',colors(j,:));
    plot(j,P.bias_change_pp(j),'o','Color',colors(j,:),'MarkerFaceColor',colors(j,:),'MarkerSize',4.5);
end
ylabel({'Bias change','(percentage points)'});ylim([-.020 .035]);
yticks(-.02:.01:.03);ax.YAxis.Exponent=0;ytickformat('%.2f');
decorate(ax,ticks,'c');

annotation(fig,'textbox',[.04 .015 .94 .125],...
    'String',{'* Clay D1   {\it h/D} < 5   Both 80-element branches failed at 10^{-9}',...
    'Diagnostic retries used a residual tolerance of 3 \times 10^{-9}'},...
    'EdgeColor','none','FontName','Times New Roman','FontSize',8.5,...
    'Interpreter','tex','HorizontalAlignment','left','VerticalAlignment','middle');
drawnow;
id='Fig36_mesh_effect_on_freezing';
savefig(fig,fullfile(out,[id '.fig']));
exportgraphics(fig,fullfile(out,[id '.pdf']),'ContentType','vector');
print(fig,fullfile(out,[id '.svg']),'-dsvg');
print(fig,fullfile(out,[id '.emf']),'-dmeta');
exportgraphics(fig,fullfile(out,[id '_QA.png']),'Resolution',180);
close(fig);

caption=['Mesh sensitivity of geometry-freezing bias at D=1.1D_0 and H equal to the original 3D load at D_0/10. '...
    'Thickness, embedment, loading height and soil inputs are fixed. '...
    '(a) Bias 100(v_F/v_U-1) for 40 and 80 embedded elements. '...
    '(b) Signed displacement changes 100(v_80/v_40-1) for updated and frozen soil diameters. '...
    '(c) The 80-element minus 40-element bias in percentage points. '...
    'The largest absolute displacement change is 0.0323% and the largest bias change is 0.0211 percentage points, '...
    'within the declared 1% and 0.5 percentage-point mesh screens. '...
    'Asterisks identify Clay D1, for which h/D=4.5455 lies below the calibration bound of 5. '...
    'Both 80-element D1 branches failed the original 10^-9 residual criterion and were retried at a diagnostic tolerance of 3 x 10^-9. '...
    'The 40-element calculations use the original 10^-9 criterion. These comparisons quantify discretization sensitivity of the specified reduced-model interventions.'];
cap=table(string(id),3,string(caption),...
    "Fig36_mesh_effect_on_freezing_plot_data.csv; Fig36_complete_response_data.csv; Fig36_retained_strict_failures.csv",...
    "Tests stability of the high-load freezing-bias claim under mesh refinement",...
    "mesh_verification/code/plot_freezing_mesh_review.m",...
    'VariableNames',{'id','panels','caption','data_sources','evidence_value','program'});
writetable(cap,fullfile(out,'figure_captions_evidence.csv'));
fid=fopen(fullfile(out,[id '_caption.txt']),'w','n','UTF-8');fprintf(fid,'%s\n',caption);fclose(fid);
summary=struct('matlab',version,'platform',computer,'group_figures',1,'panels',3,...
    'model_solves_executed',0,'existing_response_states',height(R),'paired_conditions',height(P),...
    'retained_failed_branches',height(F),'newton_tolerance_original',1e-9,...
    'D1_diagnostic_retry_tolerance',3e-9,'max_abs_displacement_change_percent',...
    max(abs([P.updated_displacement_change_percent;P.frozen_displacement_change_percent])),...
    'max_abs_bias_change_pp',max(abs(P.bias_change_pp)),...
    'nominal_mesh_screen_displacement_percent',1,'nominal_mesh_screen_bias_pp',.5);
fid=fopen(fullfile(out,'Fig36_export_record.json'),'w');fwrite(fid,jsonencode(summary,PrettyPrint=true),'char');fclose(fid);
fprintf('Exported one figure with three panels. No model solves executed.\n');
end

function decorate(ax,ticks,letter)
xlim(ax,[.6 4.4]);xticks(ax,1:4);xticklabels(ax,ticks);title(ax,'');
ax.FontSize=9;ax.XLabel.FontSize=10;ax.YLabel.FontSize=10;
text(ax,.015,.985,['(' letter ')'],'Units','normalized','VerticalAlignment','top',...
    'FontName','Times New Roman','FontSize',10.5,'FontWeight','bold','Tag','PanelLabel');
end
