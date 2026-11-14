function [y,dy]=pisa_conic(x,a)

x=x(:);s=sign(x);v=abs(x);xu=a.xu+zeros(size(x));k=a.k+zeros(size(x));n=a.n+zeros(size(x));yu=a.yu+zeros(size(x));
y=zeros(size(x));dy=y;
bil=n==0; y(bil)=min(k(bil).*v(bil),yu(bil));dy(bil)=k(bil).*(v(bil)<yu(bil)./k(bil));
lim=n==1;y(lim)=yu(lim).*min(v(lim)./xu(lim),1);dy(lim)=yu(lim)./xu(lim).*(v(lim)<xu(lim));
sat=v>=xu&~bil&~lim;y(sat)=yu(sat);
j=find(~bil&~lim&~sat);
if ~isempty(j)
 u=v(j);X=xu(j);K=k(j);N=n(j);Y=yu(j);aa=1-2*N;
 bb=2*N.*u./X-(1-N).*(1+u.*K./Y);cc=(1-N).*u.*K./Y-N.*(u./X).^2;
 disc=bb.^2-4*aa.*cc;assert(all(disc>=-1e-12));den=-bb+sqrt(max(0,disc));
 q=2*cc./den;zero=abs(den)<1e-14;q(zero)=(-bb(zero)-sqrt(max(0,disc(zero))))./(2*aa(zero));
 y(j)=Y.*q;
 bp=2*N./X-(1-N).*K./Y;cp=(1-N).*K./Y-2*N.*u./X.^2;
 dy(j)=-Y.*(bp.*q+cp)./(2*aa.*q+bb);
end
y=y.*s;

assert(all(isfinite(y))&&all(isfinite(dy))&&all(dy>=-1e-10));
end
