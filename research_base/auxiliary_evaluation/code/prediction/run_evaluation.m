function result = run_evaluation(data,outdir,cfg)

if nargin<3,cfg=struct;end
cfg=defaults(cfg);timer=tic; validate(data);
if ~isfolder(outdir),mkdir(outdir);end
assert(~isfile(fullfile(outdir,'evaluation_configuration.json')),'evaluation:ExistingRun','Choose a fresh evaluation output directory; previous runs are retained.');
[tr,ca,te,roles]=splitMasks(data.split,size(data.X,1));
group=string(data.group(:));ug=unique(group,'stable');
for k=1:numel(ug),assert(numel(unique(roles(group==ug(k))))==1,'evaluation:GroupLeakage','A geometry crosses partitions.');end
assert(nnz(tr)>0&&nnz(ca)>0&&nnz(te)>0,'evaluation:EmptyPartition','All three roles are required.');
N=size(data.X,1);J=size(data.low,2); ngcal=numel(unique(group(ca)));
if ~isfield(data,'mass'),data.mass=nan(N,1);end
data.mass=data.mass(:);
assert(numel(data.mass)==N,'evaluation:Mass','Mass needs one value per row.');
if any(isfinite(data.mass))
 assert(all(isfinite(data.mass))&&all(data.mass>0),'evaluation:Mass','Mass must be positive finite, or all omitted.');
 for k=1:numel(ug),v=data.mass(group==ug(k));assert(max(v)-min(v)<=1e-9*max(v),'evaluation:MassGroup','Mass varies within geometry.');end
end
sourceDir=fileparts(mfilename('fullpath'));
protocol=struct('config',cfg,'function_sha256',fileHash(mfilename('fullpath')),...
 'fit_sha256',fileHash(fullfile(sourceDir,'fit_response_model.m')),...
 'predict_sha256',fileHash(fullfile(sourceDir,'predict_response_model.m')),...
 'data_sha256',arrayHash(data),'created',char(datetime('now','TimeZone','Asia/Shanghai')),...
 'test_reference_used_for_tuning',false,'training_rows',nnz(tr),'calibration_rows',nnz(ca),'test_rows',nnz(te));
writeJSON(fullfile(outdir,'evaluation_configuration.json'),protocol);
writetable(table((1:N)',group,roles,'VariableNames',{'row','geometry_group','role'}),fullfile(outdir,'split_manifest.csv'));

lowmodel=fit_response_model(data.X(tr,:),data.low(tr,:),group(tr),cfg.model);
deltamodel=fit_response_model(data.X(tr,:),data.reference(tr,:)-data.low(tr,:),group(tr),cfg.model);
refmodel=fit_response_model(data.X(tr,:),data.reference(tr,:),group(tr),cfg.model);
[pl,sl,ml]=predict_response_model(lowmodel,data.X);
[pd,sd]=predict_response_model(deltamodel,data.X);
[pr,sr]=predict_response_model(refmodel,data.X);
names=["low_proxy","global_p95","direct_reference","source_corrected"];

trainGroups=unique(group(tr),'stable');groupBias=zeros(numel(trainGroups),J);
for k=1:numel(trainGroups)
 rows=group==trainGroups(k);groupBias(k,:)=mean(data.reference(rows,:)-pl(rows,:),1);
end
globalBias=mean(groupBias,1);
means=cat(3,pl,pl+globalBias,pr,pl+pd); scales=cat(3,ones(N,J),ones(N,J),sr,sl+sd);
lower=zeros(size(means));upper=lower;q=zeros(4,1);calrows=cell(4,1);
for m=1:4
 if m==1,target=data.low(ca,:);else,target=data.reference(ca,:);end
 residual=abs(target-means(ca,:,m))./scales(ca,:,m);
 [q(m),s,calgroups]=calibrate(residual,group(ca),cfg.alpha);
 lower(:,:,m)=means(:,:,m)-q(m)*scales(:,:,m);upper(:,:,m)=means(:,:,m)+q(m)*scales(:,:,m);
 calrows{m}=table(repmat(names(m),numel(s),1),calgroups,s,repmat(q(m),numel(s),1),...
  'VariableNames',{'method','geometry_group','max_nonconformity','calibration_quantile'});
end
writetable(vertcat(calrows{:}),fullfile(outdir,'calibration_scores.csv'));
writetable(lowmodel.cv_trials,fullfile(outdir,'cv_low.csv'));writetable(deltamodel.cv_trials,fullfile(outdir,'cv_delta.csv'));writetable(refmodel.cv_trials,fullfile(outdir,'cv_reference.csv'));
prediction_freeze_time=char(datetime('now','TimeZone','Asia/Shanghai'));

save(fullfile(outdir,'frozen_predictions.mat'),'names','means','scales','lower','upper','q','globalBias','tr','ca','te','group','cfg','prediction_freeze_time','-v7');
save(fullfile(outdir,'fitted_models.mat'),'lowmodel','deltamodel','refmodel','-v7');

testidx=find(te); testgroup=group(te);testug=unique(testgroup,'stable');ng=numel(testug);
truth=zeros(ng,1);point=zeros(ng,4);lo=point;hi=point;covered=false(ng,4);mass=nan(ng,1);ood=zeros(ng,1);
for k=1:ng
 rows=find(group==testug(k)); truth(k)=max(data.reference(rows,:),[],'all');mass(k)=data.mass(rows(1));ood(k)=max(ml.nearest_training_standardized_distance(rows));
 for m=1:4
  point(k,m)=max(means(rows,:,m),[],'all');lo(k,m)=max(lower(rows,:,m),[],'all');hi(k,m)=max(upper(rows,:,m),[],'all');
  covered(k,m)=all(data.reference(rows,:)>=lower(rows,:,m)&data.reference(rows,:)<=upper(rows,:,m),'all');
 end
end
groupTables=cell(4,1);rowTables=cell(4,1);summary=cell(4,1);
for m=1:4
 approved=hi(:,m)<=1;rejected=lo(:,m)>1;undecided=~(approved|rejected);bad=truth>1;
 groupTables{m}=table(repmat(names(m),ng,1),testug,mass,truth,point(:,m),lo(:,m),hi(:,m),covered(:,m),approved,rejected,undecided,ood,...
 'VariableNames',{'method','geometry_group','mass','reference_max_utilization','predicted_max_utilization','lower_max','upper_max','joint_covered','certified_feasible','certified_infeasible','unresolved','OOD_distance'});
 row=zeros(nnz(te)*J,9);at=0;
 for j=1:J
  rr=at+(1:nnz(te));row(rr,:)=[testidx,repmat(j,nnz(te),1),data.low(te,j),data.reference(te,j),means(te,j,m),scales(te,j,m),lower(te,j,m),upper(te,j,m),data.mass(te)];at=at+nnz(te);
 end
 rowTables{m}=addvars(array2table(row,'VariableNames',{'row','constraint','low','reference','prediction','scale','lower','upper','mass'}),repmat(names(m),size(row,1),1),'Before',1,'NewVariableNames','method');
 e=means(te,:,m)-data.reference(te,:);selected=find(approved);
 selectedmass=NaN;selectedvalid=NaN;
 if ~isempty(selected)&&all(isfinite(mass))
  [selectedmass,ii]=min(mass(selected));selectedvalid=double(truth(selected(ii))<=1);
 end
 summary{m}=struct('method',char(names(m)),'test_groups',ng,'reference_feasible_groups',sum(~bad),...
  'approved_groups',sum(approved),'false_approved_groups',sum(approved&bad),...
  'false_approval_given_approved',safeDivide(sum(approved&bad),sum(approved)),...
  'approval_fraction',mean(approved),'unresolved_fraction',mean(undecided),'joint_group_coverage',mean(covered(:,m)),...
  'row_response_RMSE',sqrt(mean(e.^2,'all')),'row_response_MAE',mean(abs(e),'all'),...
  'selected_certified_mass',selectedmass,'selected_reference_feasible',selectedvalid,'calibration_quantile',q(m));
end
writetable(vertcat(groupTables{:}),fullfile(outdir,'test_group_predictions.csv'));
writetable(vertcat(rowTables{:}),fullfile(outdir,'test_row_predictions.csv'));
summary=vertcat(summary{:});writetable(struct2table(summary),fullfile(outdir,'test_summary.csv'));

stream=RandStream('mt19937ar','Seed',cfg.bootstrap_seed);B=randi(stream,ng,ng,cfg.bootstrap_replicates);
F=cfg.approval_fractions(:); curves={};deltas={};ci=0;di=0;
for a=1:numel(F)
 count=max(1,floor(F(a)*ng));bootrisk=zeros(cfg.bootstrap_replicates,4);bootviol=bootrisk;risk=zeros(1,4);violation=risk;
 for m=1:4
  [risk(m),violation(m)]=rankMetrics(truth,hi(:,m),point(:,m),count);
  for b=1:cfg.bootstrap_replicates
   z=B(:,b);[bootrisk(b,m),bootviol(b,m)]=rankMetrics(truth(z),hi(z,m),point(z,m),count);
  end
  ci=ci+1;curves{ci}=struct('method',char(names(m)),'target_fraction',F(a),'approved_count',count,'approved_fraction',count/ng,...
   'false_approval_given_approved',risk(m),'false_approval_ci_low',percentile(bootrisk(:,m),.025),'false_approval_ci_high',percentile(bootrisk(:,m),.975),...
   'mean_positive_violation_approved',violation(m),'violation_ci_low',percentile(bootviol(:,m),.025),'violation_ci_high',percentile(bootviol(:,m),.975));
 end
 for baseline=[2 3]
  dr=bootrisk(:,4)-bootrisk(:,baseline);dv=bootviol(:,4)-bootviol(:,baseline);di=di+1;
  deltas{di}=struct('comparison',['source_corrected_minus_' char(names(baseline))],'target_fraction',F(a),'false_approval_difference',risk(4)-risk(baseline),...
   'difference_ci_low',percentile(dr,.025),'difference_ci_high',percentile(dr,.975),'positive_violation_difference',violation(4)-violation(baseline),...
   'violation_difference_ci_low',percentile(dv,.025),'violation_difference_ci_high',percentile(dv,.975));
 end
end
CT=struct2table(vertcat(curves{:}));DT=struct2table(vertcat(deltas{:}));writetable(CT,fullfile(outdir,'matched_approval_curves.csv'));writetable(DT,fullfile(outdir,'paired_effects.csv'));
primary=DT(abs(DT.target_fraction-.5)<1e-12&string(DT.comparison)=="source_corrected_minus_global_p95",:);
support=height(primary)==1&&primary.difference_ci_high<0&&primary.violation_difference_ci_high<=0;
run=struct('matlab_version',version,'compatibility_target','MATLAB R2024b','base_matlab_only',true,...
 'analysis_status','Executed numerical-reference validation; physical scope inherited from source model',...
 'test_reference_used_for_training_or_tuning',false,'prediction_frozen_at',prediction_freeze_time,...
 'development_groups',numel(unique(group(tr))),'calibration_groups',ngcal,'test_groups',ng,...
 'shared_reference_development_rows',nnz(tr),'shared_reference_calibration_rows',nnz(ca),'reference_test_rows_for_evaluation_only',nnz(te),...
 'low_fit_seconds',lowmodel.elapsed_seconds,'delta_fit_seconds',deltamodel.elapsed_seconds,'direct_reference_fit_seconds',refmodel.elapsed_seconds,...
 'low_fit_count',lowmodel.fit_count,'delta_fit_count',deltamodel.fit_count,'direct_reference_fit_count',refmodel.fit_count,...
 'source_corrected_development_reference_rows',nnz(tr),'global_p95_development_reference_rows',nnz(tr),...
 'global_p95_correction','Development-group-equal mean reference-minus-low-proxy residual; constant across X',...
 'reference_aware_methods_calibration_rows',nnz(ca),'low_proxy_reference_calibration_rows',0,...
 'total_elapsed_seconds',toc(timer),'bootstrap_seed',cfg.bootstrap_seed,'bootstrap_replicates',cfg.bootstrap_replicates,...
 'primary_fraction',.5,'H2_prespecified_superiority_supported',support,...
 'coverage_interpretation','Group-level empirical coverage; no guarantee over adaptive optimization or unrepresented sites',...
 'physical_validation',false);
writeJSON(fullfile(outdir,'evaluation_run.json'),run);
result=struct('summary',struct2table(summary),'curves',CT,'paired_effects',DT,'run',run,'models',struct('low',lowmodel,'delta',deltamodel,'reference',refmodel),...
 'group',testug,'truth',truth,'means',point,'lower',lo,'upper',hi,'mass',mass);
save(fullfile(outdir,'evaluation_result.mat'),'result','-v7');
end
function cfg=defaults(cfg)
d=struct('alpha',.05,'bootstrap_seed',271828,'bootstrap_replicates',500,'approval_fractions',.1:.1:1,'model',struct);
fn=fieldnames(d);for i=1:numel(fn),if ~isfield(cfg,fn{i}),cfg.(fn{i})=d.(fn{i});end,end
assert(cfg.alpha>0&&cfg.alpha<1&&cfg.bootstrap_replicates>=1&&all(cfg.approval_fractions>0&cfg.approval_fractions<=1));
end
function validate(d)
required={'X','low','reference','group','split'};for k=1:numel(required),assert(isfield(d,required{k}),'evaluation:MissingField','Missing %s.',required{k});end
assert(isnumeric(d.X)&&isnumeric(d.low)&&isnumeric(d.reference)&&isreal(d.X)&&isreal(d.low)&&isreal(d.reference));
assert(size(d.X,1)==size(d.low,1)&&isequal(size(d.low),size(d.reference))&&size(d.low,1)==numel(d.group),'evaluation:Shape','Input sizes mismatch.');
assert(all(isfinite(d.X),'all')&&all(isfinite(d.low),'all')&&all(isfinite(d.reference),'all'),'evaluation:Finite','Nonfinite rows must be resolved/recorded, not silently omitted.');
assert(all(~ismissing(string(d.group(:)))),'evaluation:Group','Missing group labels rejected.');
end
function [tr,ca,te,roles]=splitMasks(s,N)
if isstruct(s)
 roles=repmat("",N,1);fields={'dev','cal','test'};
 for i=1:3
  ix=s.(fields{i});if islogical(ix),assert(numel(ix)==N);else,assert(all(ix>=1&ix<=N&ix==floor(ix)));end
  assert(all(roles(ix)==""),'evaluation:Overlap','Partition rows overlap.');roles(ix)=fields{i};
 end
else,roles=string(s(:));end
assert(numel(roles)==N&&all(ismember(roles,["dev","cal","test"])),'evaluation:Split','Every row needs dev/cal/test role.');tr=roles=="dev";ca=roles=="cal";te=roles=="test";
end
function [q,s,groups]=calibrate(res,group,alpha)
groups=unique(group,'stable');s=zeros(numel(groups),1);
for k=1:numel(groups),s(k)=max(res(group==groups(k),:),[],'all');end
rank=ceil((numel(s)+1)*(1-alpha));ss=sort(s);if rank>numel(s),q=Inf;else,q=ss(rank);end
end
function [risk,v]=rankMetrics(truth,upper,mu,k)

[~,ix]=sortrows([upper,mu,(1:numel(truth))'],[1 2 3]);ix=ix(1:k);risk=mean(truth(ix)>1);v=mean(max(truth(ix)-1,0));
end
function y=percentile(x,p)
x=sort(x);pos=1+(numel(x)-1)*p;a=floor(pos);b=ceil(pos);y=x(a)+(pos-a)*(x(b)-x(a));
end
function y=safeDivide(a,b),if b==0,y=NaN;else,y=a/b;end,end
function writeJSON(path,value)
f=fopen(path,'w','n','UTF-8');assert(f>=0);cl=onCleanup(@()fclose(f));fwrite(f,jsonencode(value,'PrettyPrint',true),'char');
end
function h=fileHash(path)
if ~isfile(path),path=[path '.m'];end
f=fopen(path,'rb');c=onCleanup(@()fclose(f));bytes=fread(f,Inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
function h=arrayHash(value)
md=java.security.MessageDigest.getInstance('SHA-256');md.update(getByteStreamFromArray(value));h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
