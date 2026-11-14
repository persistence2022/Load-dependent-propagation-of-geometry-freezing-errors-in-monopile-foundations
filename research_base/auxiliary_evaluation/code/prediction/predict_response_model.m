function [mu,scale,meta] = predict_response_model(model,X)

assert(isnumeric(X)&&isreal(X)&&size(X,2)==size(model.Xtrain,2)&&all(isfinite(X),'all'),'evaluation:PredictionInput','Invalid prediction X.');
Xs=(X-model.xmean)./model.xstd;
d2=max(0,sum(Xs.^2,2)+sum(model.Xtrain.^2,2)'-2*(Xs*model.Xtrain'));
K=exp(-d2/(2*model.lengthscale^2*size(X,2)));
mu=model.ymean+(K*model.alpha).*model.ystd;
nominal=sqrt(max(0,1-sum((model.L\K').^2,1)'));
scale=(nominal+model.scale_floor).*model.ystd;
meta=struct('nominal_distance_scale',nominal,'nearest_training_standardized_distance',sqrt(min(d2,[],2)/size(X,2)));
end
