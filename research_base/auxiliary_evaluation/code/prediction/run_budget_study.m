function result=run_budget_study(data,c,outdir)

if ~isfolder(outdir),mkdir(outdir);end
assert(all(isfinite(data.low),'all')&&all(isfinite(data.reference),'all'));
assert(isfield(c.policy,'relative_cost'),'stage:MissingCost','Supply independently measured development query costs.');
cost=c.policy.relative_cost;
cfg=c.policy;cfg.relative_cost=cost;cfg.response_count=size(data.low,2);cfg.seed=c.seeds.policy;
stage.save_json(fullfile(outdir,'budget_protocol.json'),cfg);
rows={};allhistory={};allledger={};at=0;T0=tic;
for rep=1:cfg.repetitions
 stream=RandStream('mt19937ar','Seed',cfg.seed+rep);order=randperm(stream,size(data.X,1));initL=order(1:cfg.initial_low);initR=initL(1:cfg.initial_reference);
 for b=cfg.budget_reference_equivalents
  for method=string(cfg.methods)
   ledger=zeros(0,4); queryTimer=0;
   [selection,hist,notes]=policy(data.X,data.mass,@query,initL,initR,method,b,cfg);
   reference=data.reference(selection,:);maxU=max(reference);feasible=maxU<=1;
   at=at+1;rows{at}=struct('method',char(method),'repeat',rep,'additional_budget_R_equivalent',b,'selected_geometry',selection,...
    'mass_kg',data.mass(selection),'reference_max_utilization',maxU,'reference_feasible',feasible,'positive_violation',max(maxU-1,0),...
    'charged_query_cost_R_equivalent',sum(ledger(:,3)),'initial_cost_R_equivalent',numel(initL)*cost(1)+numel(initR),...
    'online_cost_R_equivalent',notes.spent,'low_queries',sum(ledger(:,2)==1),'reference_queries',sum(ledger(:,2)==2),...
    'final_verification_cost_R_equivalent',1,'policy_fit_seconds',notes.fit_seconds,'lookup_wall_seconds',queryTimer,...
    'recommendation_type',notes.recommendation_type);
   H=array2table(hist,'VariableNames',{'iteration','spent_R_equivalent','candidate','action','selected_mass_kg','estimated_max_utilization','known_reference_feasible_best_mass_kg'});
   H=addvars(H,repmat(method,height(H),1),repmat(rep,height(H),1),repmat(b,height(H),1),'Before',1,'NewVariableNames',{'method','repeat','budget'});allhistory{at}=H;
   L=array2table(ledger,'VariableNames',{'geometry_id','level','charged_cost_R_equivalent','lookup_seconds'});
   L=addvars(L,repmat(method,height(L),1),repmat(rep,height(L),1),repmat(b,height(L),1),'Before',1,'NewVariableNames',{'method','repeat','budget'});allledger{at}=L;
  end
 end
 fprintf('Budget study %d/%d initializations complete\n',rep,cfg.repetitions);
end
T=struct2table(vertcat(rows{:}));valid=max(data.reference,[],2)<=1;
if any(valid),best=min(data.mass(valid));else,best=NaN;end
T.reference_candidate_optimum_kg=repmat(best,height(T),1);T.feasible_mass_regret_fraction=(T.mass_kg-best)/best;T.feasible_mass_regret_fraction(~T.reference_feasible)=NaN;
writetable(T,fullfile(outdir,'policy_results.csv'));writetable(vertcat(allhistory{:}),fullfile(outdir,'policy_history.csv'));writetable(vertcat(allledger{:}),fullfile(outdir,'query_ledger.csv'));
groups=unique(T(:,{'method','additional_budget_R_equivalent'}),'rows','stable');summary=cell(height(groups),1);
for i=1:height(groups)

 methodName=string(groups.method(i,:));
 r=T(string(T.method)==methodName&T.additional_budget_R_equivalent==groups.additional_budget_R_equivalent(i),:);v=r.reference_feasible;
 summary{i}=struct('method',char(methodName),'budget',groups.additional_budget_R_equivalent(i),'runs',height(r),'feasible_runs',sum(v),...
 'feasible_fraction',mean(v),'mean_positive_violation',mean(r.positive_violation),'mean_conditional_mass_kg',mean(r.mass_kg(v)),...
 'mean_conditional_regret',mean(r.feasible_mass_regret_fraction(v)),'mean_total_query_cost',mean(r.charged_query_cost_R_equivalent+1),'mean_fit_seconds',mean(r.policy_fit_seconds));
end
S=struct2table(vertcat(summary{:}));writetable(S,fullfile(outdir,'policy_summary.csv'));
record=struct('actual_matlab',version,'candidate_geometries',size(data.X,1),'initializations',cfg.repetitions,'simulations',height(T),...
 'candidate_reference_feasible',sum(valid),'candidate_reference_best_mass_kg',best,'low_to_reference_measured_cost_ratio',cost(1),...
 'execution_seconds',toc(T0),'independence','Repeated algorithm initialization on one numerical candidate domain; not independent sites or physical trials',...
 'offline_table_generation_seconds',sum(data.elapsed,'all'),'unused_models_in_offline_generation','updated_linear and frozen_nonlinear retained for H1; not charged as policy queries',...
 'information_control','policy receives X,mass,query handle only; queried labels and costs logged; final reference truth outside policy',...
 'scale_interpretation','Heuristic nominal kernel scales; not calibrated safety probabilities','physical_validation',false);
stage.save_json(fullfile(outdir,'policy_run.json'),record);result=struct('table',T,'summary',S,'run',record);save(fullfile(outdir,'policy_result.mat'),'result','-v7');
 function y=query(index,level)
  assert(isscalar(index)&&index>=1&&index<=size(data.X,1)&&ismember(level,[1 2]));tm=tic;
  if level==1,y=data.low(index,:);else,y=data.reference(index,:);end
  dt=toc(tm);queryTimer=queryTimer+dt;ledger(end+1,:)=[index level cost(level) dt];
 end
end

function [selected,history,note]=policy(X,mass,query,initialL,initialR,method,budget,c)
n=size(X,1);J=c.response_count;YL=nan(n,J);YR=YL;knownL=false(n,1);knownR=knownL;
for i=initialL,YL(i,:)=query(i,1);knownL(i)=true;end
for i=initialR,YR(i,:)=query(i,2);knownR(i)=true;end
fitcfg=struct('tune',false,'lengthscale',c.fixed_kernel_length,'ridge',c.ridge,'scale_floor',.02);
spent=0;fitseconds=0;iteration=0;history=zeros(0,7);type='predicted feasible with heuristic margin';
while true
 tm=tic;
 if method=="ReferenceOnly"
  refModel=fit_response_model(X(knownR,:),YR(knownR,:),find(knownR),fitcfg);[mu,sigma]=predict_response_model(refModel,X);sL=zeros(n,J);sD=sigma;
 else
  lm=fit_response_model(X(knownL,:),YL(knownL,:),find(knownL),fitcfg);dm=fit_response_model(X(knownR,:),YR(knownR,:)-YL(knownR,:),find(knownR),fitcfg);
  [pl,sL]=predict_response_model(lm,X);[pd,sD]=predict_response_model(dm,X);mu=pl+pd;sigma=sL+sD;
 end
 fitseconds=fitseconds+toc(tm);mu(knownR,:)=YR(knownR,:);sigma(knownR,:)=0;
 U=max(mu+c.scale_multiplier*sigma,[],2);G=max(mu,[],2);s=max(sigma,[],2);feasible=U<=1;
 certifiedKnown=knownR&max(YR,[],2)<=1;inc=Inf;if any(certifiedKnown),inc=min(mass(certifiedKnown));end
 if any(feasible),idx=find(feasible);[~,j]=min(mass(idx));selected=idx(j);type='predicted feasible with heuristic margin';
 elseif any(certifiedKnown),idx=find(certifiedKnown);[~,j]=min(mass(idx));selected=idx(j);type='queried reference feasible';
 else,[~,selected]=min(G);type='least predicted utilization; feasibility unresolved';end
 if knownR(selected)&&max(YR(selected,:))<=1,type='queried reference feasible';end
 if method=="Frozen"||spent>=budget-1e-12||iteration>=2*n,break;end
 if isfinite(inc),gain=max(inc-mass,0)/inc;else,gain=ones(n,1);end

 p=0.5*erfc((G-1)./(sqrt(2)*max(s,1e-4)));boundary=exp(-0.5*((G-1)./max(s,1e-3)).^2);
 weight=.65*gain.*p+.30*boundary+.05*s/max(max(s),1e-9);
 scoreL=weight.*max(sL,[],2)/c.relative_cost(1);scoreR=weight.*max(sD,[],2);
 scoreL(knownL)=-Inf;scoreR(knownR)=-Inf;
 if budget-spent<c.relative_cost(1)-1e-12,scoreL(:)=-Inf;end
 highCost=ones(n,1);if method~="ReferenceOnly",highCost=highCost+c.relative_cost(1)*~knownL;end
 scoreR=scoreR./highCost;scoreR(highCost>budget-spent+1e-12)=-Inf;
 if method=="LowOnly",scoreR(:)=-Inf;elseif method=="ReferenceOnly",scoreL(:)=-Inf;
 elseif method=="Periodic"
  if mod(iteration,3)==2&&any(isfinite(scoreR)),scoreL(:)=-Inf;elseif any(isfinite(scoreL)),scoreR(:)=-Inf;end
 end
 [value,choice]=max([scoreL;scoreR]);if ~isfinite(value),break;end
 level=1+(choice>n);index=mod(choice-1,n)+1;
 if level==2&&method~="ReferenceOnly"&&~knownL(index),YL(index,:)=query(index,1);knownL(index)=true;spent=spent+c.relative_cost(1);end
 if level==1,YL(index,:)=query(index,1);knownL(index)=true;spent=spent+c.relative_cost(1);else,YR(index,:)=query(index,2);knownR(index)=true;spent=spent+1;end
 iteration=iteration+1;history(end+1,:)=[iteration spent index level mass(selected) G(selected) inc];
end
assert(spent<=budget+1e-9);note=struct('spent',spent,'fit_seconds',fitseconds,'recommendation_type',type);
end
