function report=run_supplement(mode)
if nargin==0,mode='checks';end
mode=validatestring(mode,{'checks','smoke','audit','figures','formal','generate_auxiliary'});
root=fileparts(mfilename('fullpath'));oldPath=path;oldDirectory=pwd;
restore=onCleanup(@()restoreSession(oldPath,oldDirectory));
base=fullfile(root,'research_base');research=fullfile(base,'research');
addpath(research,fullfile(research,'code'),fullfile(research,'pisa_model'),'-begin');
out=fullfile(root,'reproduction_runs',[mode '_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);mkdir(out);
report=struct('mode',mode,'started',char(datetime('now')),'matlab',version,'release',version('-release'),'platform',computer,'status','RUNNING');
timer=tic;
try
    switch mode
        case {'checks','smoke'}
            ec=fullfile(base,'auxiliary_evaluation','code');addpath(ec,'-begin');verify_portable_inputs;
            report.execution=run_research(mode);
        case 'audit'
            ec=fullfile(base,'auxiliary_evaluation','code');addpath(ec,'-begin');verify_portable_inputs;run_evaluation_audit;
            report.execution=jsondecode(fileread(fullfile(base,'auxiliary_evaluation','current_audit.json')));
        case 'figures'
            report.execution=run_research('figures');
            plots=fullfile(base,'auxiliary_evaluation','plotting');
            addpath(fullfile(plots,'auxiliary_plots'),fullfile(plots,'evaluation_plots'),'-begin');
            plot_auxiliary_results(fullfile(plots,'reference_results'));
            plot_evaluation_results(fullfile(plots,'reference_diagnostics'));
            addpath(fullfile(base,'mesh_verification','code'),'-begin');plot_freezing_mesh_review;
        case 'formal'
            report.execution=run_research('all');
        case 'generate_auxiliary'
            aux=fullfile(root,'auxiliary_generation');addpath(aux,fullfile(aux,'physics'),'-begin');cd(aux);
            c=analysis_config();fresh=fullfile(aux,'runs',['generated_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);mkdir(fresh);
            stage.save_json(fullfile(fresh,'configuration.json'),c);
            [data,info]=stage.generate_data(c,'formal',fresh);assert(info.invalid_model_rows==0);
            [pool,poolinfo]=stage.generate_data(c,'policy',fresh);assert(poolinfo.invalid_model_rows==0);
            evaluation=fullfile(base,'auxiliary_evaluation','code','prediction');addpath(evaluation,'-begin');
            opts=struct('approval_fractions',unique([.1:.1:1,.25,.75]),'bootstrap_seed',271828);
            run_evaluation(data,fullfile(fresh,'evaluation'),opts);
            dev=string(data.split)=="dev";costs=[median(data.elapsed(dev,1)),median(data.elapsed(dev,4))];
            c.policy.relative_cost=costs/costs(2);c.policy.reference_solve_seconds=costs(2);
            stage.save_json(fullfile(fresh,'development_query_cost.json'),struct('cost_seconds',costs,'relative_cost',c.policy.relative_cost,'source','Development paired-load solves'));
            run_budget_study(pool,c,fullfile(fresh,'policy'));
            report.execution=struct('run',fresh,'geometry_generation',info,'pool_generation',poolinfo,'scope','Fresh numerical responses and timing-dependent policy calculations');
    end
    report.status='PASS';report.elapsed_seconds=toc(timer);saveJSON(fullfile(out,'runtime_report.json'),report);
    fprintf('%s PASS in %.4f seconds\n',mode,report.elapsed_seconds);
catch err
    report.status='FAILED';report.elapsed_seconds=toc(timer);report.error=getReport(err,'extended','hyperlinks','off');saveJSON(fullfile(out,'runtime_report.json'),report);rethrow(err);
end
end
function restoreSession(p,d),path(p);cd(d);end
function saveJSON(p,v)
fid=fopen(p,'w','n','UTF-8');assert(fid>=0);cl=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(v,'PrettyPrint',true));
end
