function X=designs(n,bounds,seed)
stream=RandStream('mt19937ar','Seed',seed);d=size(bounds,1);X=zeros(n,d);
for j=1:d
 z=(randperm(stream,n)'-rand(stream,n,1))/n;X(:,j)=bounds(j,1)+z*(bounds(j,2)-bounds(j,1));
end
end
