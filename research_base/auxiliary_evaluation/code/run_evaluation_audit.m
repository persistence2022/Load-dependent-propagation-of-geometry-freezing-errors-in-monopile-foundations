function run_evaluation_audit

here=fileparts(mfilename('fullpath'));review=fileparts(here);

source=fullfile(review,'inputs','formal');

assert(isfolder(source),'review:MissingInputPackage', ...
 'The complete auxiliary_evaluation/inputs folder must accompany this code.');

addpath(fullfile(here,'prediction'));
out=fullfile(review,'runs',['audit_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);mkdir(out);
diary(fullfile(out,'matlab_run.log'));cl=onCleanup(@()diary('off'));
fprintf('Actual MATLAB %s %s\n',version,computer);
files={fullfile(source,'data_formal.mat'),fullfile(source,'data_policy.mat'),fullfile(source,'evaluation','frozen_predictions.mat'),fullfile(source,'evaluation','fitted_models.mat'),fullfile(source,'policy','policy_results.csv'),fullfile(source,'policy','query_ledger.csv'),fullfile(source,'configuration.json'),fullfile(source,'development_query_cost.json')};
before=cellfun(@stage.hash_file,files,'UniformOutput',false);
F=load(files{1});P=load(files{2});A=load(files{3});M=load(files{4});d=F.data;pd=P.data;
c=jsondecode(fileread(files{7}));cost=jsondecode(fileread(files{8}));

c.policy.methods=reshape(c.policy.methods,1,[]);
c.policy.budget_reference_equivalents=reshape(c.policy.budget_reference_equivalents,1,[]);
T=readtable(files{5},'TextType','string');Q=readtable(files{6},'TextType','string');
checks={};details=struct;tr=string(d.split)=="dev";ca=string(d.split)=="cal";te=string(d.split)=="test";
check('Partition sizes',isequal([sum(tr),sum(ca),sum(te)],[80 60 120]),'One row per geometry and max over two prescribed load cases');
check('Unique physical geometry keys',height(unique(array2table(d.X),'rows'))==260,'Full D t L key is unique');
check('Group partitions disjoint',isempty(intersect(d.group(tr),d.group(ca)))&&isempty(intersect(d.group(tr),d.group(te)))&&isempty(intersect(d.group(ca),d.group(te))),'Group identity checked before model fits');
check('Policy geometries separate',isempty(intersect(d.X,pd.X,'rows')),'The separate finite policy pool does not duplicate H2 design coordinates');
check('No failed response omitted',all(isfinite(d.low),'all')&&all(isfinite(d.reference),'all')&&all(strlength(string(d.failures))==0,'all'),'No nonfinite or failed labels in the formal evaluation');
check('Utilization and target mapping',max(abs(d.low-squeeze(d.utilization_all(:,1,:))),[],'all')<1e-14&&max(abs(d.reference-squeeze(d.utilization_all(:,4,:))),[],'all')<1e-14,'Low is frozen linear and reference is updated nonlinear');
check('Response normalization',max(abs(d.utilization_all-d.responses./reshape(d.limits,260,1,3)),[],'all')<1e-14,'Constraint limits applied before response comparison');
check('Training preprocessing',max(abs(M.lowmodel.xmean-mean(d.X(tr,:),1)))<1e-13&&max(abs(M.lowmodel.xstd-std(d.X(tr,:),0,1)))<1e-13,'Development coordinates only');
check('Equal CV budgets',M.lowmodel.fit_count==31&&M.deltamodel.fit_count==31&&M.refmodel.fit_count==31,'Six parameter pairs times five folds plus one final fit per response model');
models={M.lowmodel,M.deltamodel,M.refmodel};targets={d.low(tr,:),d.reference(tr,:)-d.low(tr,:),d.reference(tr,:)};
cvdiff=zeros(3,1);
for k=1:3
 m=models{k};Y=targets{k};loss=zeros(height(m.cv_trials),1);
 for j=1:height(m.cv_trials)
  vals=zeros(5,1);
  for fold=1:5
   it=m.training_fold~=fold;iv=~it;X=d.X(tr,:);xt=X(it,:);yt=Y(it,:);
   xm=mean(xt);xs=std(xt);xs(xs<1e-12)=1;ym=mean(yt);ys=std(yt);ys(ys<1e-12)=1;
   xt=(xt-xm)./xs;xv=(X(iv,:)-xm)./xs;ell=m.cv_trials.lengthscale(j);ridge=m.cv_trials.ridge(j);
   K=exp(-dist2(xt,xt)/(2*ell^2*size(X,2)));Ks=exp(-dist2(xv,xt)/(2*ell^2*size(X,2)));
   mu=ym+(Ks*((K+ridge*eye(size(K)))\((yt-ym)./ys))).*ys;
   vals(fold)=mean(((mu-Y(iv,:))./ys).^2,'all');
  end
  loss(j)=mean(vals);
 end
 cvdiff(k)=max(abs(loss-m.cv_trials.group_equal_scaled_MSE));
end
check('Independent fold-local CV reconstruction',max(cvdiff)<1e-7,'Independent matrix solve with train-fold means and scales matches saved CV objectives');details.CV_max_absolute_difference=cvdiff;
ec=A.cfg;run_evaluation(d,fullfile(out,'h2_replay'),ec);B=load(fullfile(out,'h2_replay','frozen_predictions.mat'));
check('Frozen H2 prediction replay',max(abs(B.means-A.means),[],'all')<1e-11&&max(abs(B.upper-A.upper),[],'all')<1e-10,'Unchanged frozen source code and original observations');
e=d;e.reference(te,:)=e.reference(te,:)+10;run_evaluation(e,fullfile(out,'h2_test_label_poison'),ec);C=load(fullfile(out,'h2_test_label_poison','frozen_predictions.mat'));
check('Test response isolation',isequal(B.means,C.means)&&isequal(B.scales,C.scales)&&isequal(B.q,C.q),'Large deterministic change to every test reference leaves fits and calibrated predictions identical');
e=d;e.reference(ca,:)=e.reference(ca,:)+7;run_evaluation(e,fullfile(out,'h2_calibration_label_poison'),ec);C=load(fullfile(out,'h2_calibration_label_poison','frozen_predictions.mat'));
check('Calibration isolated from fitting',isequal(B.means,C.means)&&isequal(B.scales,C.scales)&&any(B.q~=C.q),'Calibration changes intervals but not fitted response models');
scores=readtable(fullfile(source,'evaluation','calibration_scores.csv'),'TextType','string');
for k=1:4
 s=sort(scores.max_nonconformity(scores.method==string(A.names(k))));
 check(['Finite calibration quantile ' char(A.names(k))],numel(s)==60&&abs(s(ceil(61*.95))-A.q(k))<1e-12,'Rank 58 of 60 group-max residuals');
end
old=readtable(fullfile(source,'evaluation','paired_effects.csv'),'TextType','string');new=readtable(fullfile(out,'h2_replay','paired_effects.csv'),'TextType','string');
check('Paired group bootstrap replay',max(abs(old{:,3:end}-new{:,3:end}),[],'all')<1e-13,'500 identical resamples over 120 geometry groups');
check('Additive mean null control',max(abs(A.means(:,:,3)-A.means(:,:,4)),[],'all')<1e-10,'Same selected kernel and training coordinates make the two predictive means algebraically equal');
c.policy.relative_cost=cost.relative_cost(:)';c.policy.reference_solve_seconds=cost.cost_seconds(2);
run_budget_study(pd,c,fullfile(out,'h3_replay'));RT=readtable(fullfile(out,'h3_replay','policy_results.csv'),'TextType','string');RQ=readtable(fullfile(out,'h3_replay','query_ledger.csv'),'TextType','string');
check('H3 policy selections replay',isequal(T.selected_geometry,RT.selected_geometry)&&max(abs(T.mass_kg-RT.mass_kg))<1e-7,'All 120 strategy runs reproduce selections on the same finite pool');
check('H3 query ledger replay',isequal(Q(:,1:6),RQ(:,1:6)),'Query sequence and charges reproduced excluding wall timing');
check('Development cost source',max(abs([median(d.elapsed(tr,1)),median(d.elapsed(tr,4))]-cost.cost_seconds(:)'))<1e-14,'Cost ratio is independent development timing rather than policy lookup timing');
ledgerOK=true;extraL=0;finalKnown=0;
for k=1:height(T)
 z=T(k,:);q=Q(Q.method==z.method&Q.repeat==z.repeat&Q.budget==z.additional_budget_R_equivalent,:);
 initCost=12*c.policy.relative_cost(1)+4;
 ledgerOK=ledgerOK&&abs(sum(q.charged_cost_R_equivalent)-z.charged_query_cost_R_equivalent)<1e-10&&abs(z.initial_cost_R_equivalent-initCost)<1e-10&&abs(z.charged_query_cost_R_equivalent-initCost-z.online_cost_R_equivalent)<1e-10&&z.online_cost_R_equivalent<=z.additional_budget_R_equivalent+1e-10&&height(unique(q(:,{'geometry_id','level'}),'rows'))==height(q)&&z.final_verification_cost_R_equivalent==1;
 knownL=false(128,1);knownR=knownL;
 for j=1:height(q)
  id=q.geometry_id(j);if q.level(j)==1,knownL(id)=true;if j>16&&j<height(q)&&q.level(j+1)==2&&q.geometry_id(j+1)==id,extraL=extraL+1;end
  else,if z.method~="ReferenceOnly",ledgerOK=ledgerOK&&knownL(id);end;knownR(id)=true;end
 end
 finalKnown=finalKnown+knownR(z.selected_geometry);
end
check('H3 budgets pairing and logical final fees',ledgerOK,'Initialization and required low prequeries charged and all action budgets respected');
details.H3_immediate_low_before_reference=extraL;details.H3_final_reference_already_queried=finalKnown;

cc=c;cc.policy.repetitions=1;cc.policy.budget_reference_equivalents=16;cc.policy.methods={'DecisionAdaptive'};
qq=Q(Q.method=="DecisionAdaptive"&Q.repeat==1&Q.budget==16,:);poison=pd;
knownL=ismember((1:128)',qq.geometry_id(qq.level==1));knownR=ismember((1:128)',qq.geometry_id(qq.level==2));
poison.low(~knownL,:)=poison.low(~knownL,:)+37;poison.reference(~knownR,:)=poison.reference(~knownR,:)+53;
run_budget_study(poison,cc,fullfile(out,'h3_unqueried_label_poison'));PQ=readtable(fullfile(out,'h3_unqueried_label_poison','query_ledger.csv'),'TextType','string');PT=readtable(fullfile(out,'h3_unqueried_label_poison','policy_results.csv'),'TextType','string');
expected=T(T.method=="DecisionAdaptive"&T.repeat==1&T.additional_budget_R_equivalent==16,:);
check('Adaptive unqueried response isolation',isequal(qq(:,1:6),PQ(:,1:6))&&PT.selected_geometry==expected.selected_geometry,'Unqueried labels changed while query order and recommendation stayed fixed');
after=cellfun(@stage.hash_file,files,'UniformOutput',false);check('Reference inputs preserved',isequal(before,after),'Reference scientific arrays preserved');
writetable(cell2table(checks,'VariableNames',{'check','passed','summary'}),fullfile(out,'audit_checks.csv'));
writetable(table(string(files)',string(before)',string(after)','VariableNames',{'source','sha_before','sha_after'}),fullfile(out,'input_integrity.csv'));
record=struct('matlab',version,'platform',computer,'source_run',source, ...
 'input_package','Local reference arrays in inputs/formal', ...
 'execution_scope','Internal replay and post-evaluation adversarial diagnostics only','checks',size(checks,1),'passed',sum(cell2mat(checks(:,2))),'details',details,'original_files_unchanged',isequal(before,after));
stage.save_json(fullfile(out,'audit_summary.json'),record);stage.save_json(fullfile(review,'current_audit.json'),struct('run',out,'summary',record));
disp(jsonencode(record,'PrettyPrint',true));assert(all(cell2mat(checks(:,2))),'review:FailedCheck','See audit_checks.csv');
 function check(name,passed,evidence)
  checks(end+1,:)={name,logical(passed),evidence};fprintf('%s %s\n',string(passed),name);
 end
end
function v=dist2(a,b),v=max(0,sum(a.^2,2)+sum(b.^2,2)'-2*a*b');end
