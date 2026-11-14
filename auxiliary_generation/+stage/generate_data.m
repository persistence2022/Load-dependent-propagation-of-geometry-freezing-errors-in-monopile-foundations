function [data,info]=generate_data(c,phase,outdir)
root=fileparts(fileparts(mfilename('fullpath')));phys=fullfile(root,'physics');addpath(phys);
pc=static_config();pc.max_element_m=c.max_element_m;
sources={'static_config.m','solve_static.m','soil_py.m'};fingerprint='';
for k=1:numel(sources),fingerprint=[fingerprint sources{k} ':' stage.hash_file(fullfile(phys,sources{k})) ';'];end
physicsHash=stage.hash_bytes(unicode2native(fingerprint,'UTF-8'));cacheDir=fullfile(root,'cache','static');if ~isfolder(cacheDir),mkdir(cacheDir);end
if strcmp(phase,'formal'),roles={'development','calibration','test'};elseif strcmp(phase,'mechanism'),roles={};elseif strcmp(phase,'policy'),roles={'policy'};else,roles={phase};end
X=[];split=strings(0,1);
for r=1:numel(roles)
 role=roles{r};XX=stage.designs(c.sizes.(role),c.bounds,c.seeds.(role));X=[X;XX];split=[split;repmat(string(role),size(XX,1),1)];
end
if strcmp(phase,'mechanism')
 [d,t,l]=ndgrid(c.mechanism_grid.D_m,c.mechanism_grid.t_m,c.mechanism_grid.L_m);X=[d(:),t(:),l(:)];split=repmat("mechanism",size(X,1),1);
end
n=size(X,1);raw=nan(n,numel(c.model_modes),3);elapsed=nan(n,numel(c.model_modes));iterations=elapsed;cached=false(n,numel(c.model_modes));failures=strings(n,numel(c.model_modes));
total=tic;
for i=1:n
 for a=1:numel(c.model_modes)
  mode=c.model_modes{a};input=struct('X',X(i,:),'mode',mode,'loads',c.loads_N_Nm,'physics_configuration',pc,'physics_hash',physicsHash);
  key=stage.hash_bytes(unicode2native(jsonencode(input),'UTF-8'));cache=fullfile(cacheDir,[key '.mat']);meta=fullfile(cacheDir,[key '.json']);hit=false;
  if isfile(cache)&&isfile(meta)
   m=jsondecode(fileread(meta));if strcmp(m.key,key)&&strcmp(m.output_sha256,stage.hash_file(cache)),saved=load(cache,'result');result=saved.result;hit=true;end
  end
  if ~hit
   rt=tic;vals=nan(size(c.loads_N_Nm,1),3);it=0;status=true;msgs=strings(0);records=cell(size(c.loads_N_Nm,1),1);
   for j=1:size(c.loads_N_Nm,1)
    pp=pc;pp.H_N=c.loads_N_Nm(j,1);pp.M_Nm=c.loads_N_Nm(j,2);
    try
     q=solve_static(pp,X(i,1),X(i,2),X(i,3),mode);records{j}=q;
     if ~q.converged || any(~isfinite([q.rotation_rad q.displacement_m q.max_stress_Pa])),status=false;msgs(end+1)="nonconverged or invalid";else,vals(j,:)=[abs(q.rotation_rad),abs(q.displacement_m),q.max_stress_Pa];end
     it=it+q.iterations;
    catch err,status=false;msgs(end+1)=string(err.identifier)+": "+string(err.message);records{j}=struct('error',char(msgs(end)));end
   end
   result=struct('key',key,'input',input,'responses',max(vals,[],1),'case_responses',vals,'iterations',it,'seconds',toc(rt),'success',status,'failure_messages',{cellstr(msgs)},'case_records',{records});
   temp=[cache '.tmp.mat'];save(temp,'result','-v7');movefile(temp,cache,'f');stage.save_json(meta,struct('key',key,'output_sha256',stage.hash_file(cache)));
  end
  if result.success,raw(i,a,:)=result.responses;else,failures(i,a)=strjoin(string(result.failure_messages),'; ');end
  elapsed(i,a)=result.seconds;iterations(i,a)=result.iterations;cached(i,a)=hit;
 end
 if mod(i,20)==0||i==n,fprintf('Static %s %d/%d completed, wall %.2f s\n',phase,i,n,toc(total));end
end
limits=[c.rotation_limit_rad*ones(n,1),X(:,3)*c.displacement_length_ratio,c.stress_limit_Pa*ones(n,1)];U=raw./reshape(limits,n,1,3);
data=struct('X',X,'low',squeeze(U(:,1,:)),'reference',squeeze(U(:,4,:)),'responses',raw,'utilization_all',U,...
 'limits',limits,'mass',c.density_kg_m3*pi.*(X(:,1).^2-(X(:,1)-2*X(:,2)).^2)/4.*X(:,3),...
 'split',split,'group',(1:n)','elapsed',elapsed,'cached',cached,'failures',failures,'model_modes',{c.model_modes});
data.split(data.split=="development")="dev";data.split(data.split=="calibration")="cal";
info=struct('phase',phase,'physics_hash',physicsHash,'geometries',n,'model_pairs_per_geometry',4,'load_cases',size(c.loads_N_Nm,1),...
 'cache_hits',nnz(cached),'cache_misses',numel(cached)-nnz(cached),'invalid_model_rows',nnz(strlength(failures)>0),'actual_wall_seconds',toc(total),...
 'median_model_seconds',median(elapsed,1),'total_original_model_seconds',sum(elapsed,'all'),'scope',c.scope);
save(fullfile(outdir,['data_' phase '.mat']),'data','info','-v7');stage.save_json(fullfile(outdir,['data_' phase '_record.json']),info);
T=table(data.group,split,X(:,1),X(:,2),X(:,3),data.mass,'VariableNames',{'geometry_id','split','D_m','t_m','L_m','mass_kg'});
for a=1:4,for j=1:3,T.([c.model_modes{a} '_q' num2str(j)])=raw(:,a,j);T.([c.model_modes{a} '_util' num2str(j)])=U(:,a,j);end,end
writetable(T,fullfile(outdir,['data_' phase '.csv']));
F=table(repelem((1:n)',4),repmat(string(c.model_modes(:)),n,1),reshape(failures.',[],1),'VariableNames',{'geometry_id','model','failure'});writetable(F,fullfile(outdir,['failures_' phase '.csv']));
end
