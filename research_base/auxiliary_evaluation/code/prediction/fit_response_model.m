function model = fit_response_model(X,Y,group,cfg)

if nargin<4,cfg=struct;end
cfg=defaults(cfg);check(X,Y,group);
group=string(group(:)); groups=unique(group,'stable');
if numel(groups)<3,error('evaluation:FewGroups','At least 3 development groups are required.');end
timer=tic; trials=zeros(0,3);
if cfg.tune
 stream=RandStream('mt19937ar','Seed',cfg.seed); order=randperm(stream,numel(groups));
 fold=zeros(numel(groups),1);fold(order)=mod((0:numel(groups)-1)',min(cfg.folds,numel(groups)))+1;
 [~,gi]=ismember(group,groups); rf=fold(gi); nfold=max(fold); counts=0;
 for ell=cfg.lengthscales
  for ridge=cfg.ridges
   foldloss=zeros(nfold,1);
   for k=1:nfold
    tr=rf~=k;va=~tr; m=fitFixed(X(tr,:),Y(tr,:),ell,ridge,cfg.scale_floor);
    mu=predict_response_model(m,X(va,:));
    err=mean(((mu-Y(va,:))./m.ystd).^2,2);
    vg=group(va);ug=unique(vg); gl=zeros(numel(ug),1);
    for j=1:numel(ug),gl(j)=mean(err(vg==ug(j)));end
    foldloss(k)=mean(gl);counts=counts+1;
   end
   trials(end+1,:)=[ell,ridge,mean(foldloss)];
  end
 end
 [~,best]=min(trials(:,3));ell=trials(best,1);ridge=trials(best,2);
else
 ell=cfg.lengthscale;ridge=cfg.ridge;counts=0;rf=zeros(size(group));
end
model=fitFixed(X,Y,ell,ridge,cfg.scale_floor);
model.cv_trials=array2table(trials,'VariableNames',{'lengthscale','ridge','group_equal_scaled_MSE'});
model.training_group=group;model.training_fold=rf;model.train_rows=size(X,1);model.train_groups=numel(groups);
model.fit_count=counts+1;model.elapsed_seconds=toc(timer);model.config=cfg;
end
function m=fitFixed(X,Y,ell,ridge,floorValue)
m.xmean=mean(X,1);m.xstd=std(X,0,1);m.xstd(m.xstd<1e-12)=1;
m.ystd=std(Y,0,1);m.ystd(m.ystd<1e-12)=1;m.ymean=mean(Y,1);
m.Xtrain=(X-m.xmean)./m.xstd;Yn=(Y-m.ymean)./m.ystd;
n=size(X,1);p=size(X,2);d2=max(0,sum(m.Xtrain.^2,2)+sum(m.Xtrain.^2,2)'-2*(m.Xtrain*m.Xtrain'));
K=exp(-d2/(2*ell^2*p));m.L=chol(K+ridge*eye(n),'lower');m.alpha=m.L'\(m.L\Yn);
m.lengthscale=ell;m.ridge=ridge;m.scale_floor=floorValue;m.kind='Gaussian kernel ridge; nominal interpolation scale only';
end
function cfg=defaults(cfg)
d=struct('tune',true,'seed',314159,'folds',5,'lengthscales',[.6 1.2 2.4],...
 'ridges',[1e-6 1e-3],'lengthscale',1.2,'ridge',1e-6,'scale_floor',.02);
names=fieldnames(d);for i=1:numel(names),if ~isfield(cfg,names{i}),cfg.(names{i})=d.(names{i});end,end
assert(all(isfinite([cfg.lengthscales,cfg.ridges,cfg.lengthscale,cfg.ridge,cfg.scale_floor]))&&all([cfg.lengthscales,cfg.ridges,cfg.lengthscale,cfg.ridge,cfg.scale_floor]>0),'evaluation:Hyperparameters','Positive finite kernel settings required.');
assert(cfg.folds>=2&&cfg.folds==floor(cfg.folds),'evaluation:Folds','At least two folds required.');
end
function check(X,Y,group)
assert(isnumeric(X)&&isnumeric(Y)&&ismatrix(X)&&ismatrix(Y)&&isreal(X)&&isreal(Y),'evaluation:Type','Real numeric X/Y required.');
assert(size(X,1)==size(Y,1)&&numel(group)==size(X,1)&&size(X,2)>0&&size(Y,2)>0,'evaluation:Shape','Inconsistent X/Y/group sizes.');
assert(all(isfinite(X),'all')&&all(isfinite(Y),'all'),'evaluation:Finite','Nonfinite model data rejected.');
end
