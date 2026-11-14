function a=pisa_parameters(model,z,D,L,Dr)

z=z(:);r=z/D;s=z/L;ld=L/D;
if lower(string(model))=="cowden"
 a.p=pack(241.4+0*z,10.60-1.650*r,.9390-.03345*r,10.70-7.101*exp(-.3085*r));
 km=1.420-.09643*r;mu=.2899-.04775*r;a.m=pack(mu./km,km,0*z,mu);
 a.H=pack(235.7,2.717-.3575*ld,.8793-.03150*ld,.4038+.04812*ld);
 a.M=pack(173.1,.2146-.002132*ld,1.079-.1087*ld,.8192-.08588*ld);
else
 Dr=Dr(:);if isscalar(Dr),Dr=Dr+zeros(size(z));end
 a.p=pack(146.1-92.11*Dr,8.731-.6982*Dr-.9178*r,.917+.06193*Dr,.3667+25.89*Dr+(.3375-8.900*Dr).*s);
 mu=.2605+(-.1989+.2019*Dr).*s;a.m=pack(mu/17,17+0*z,0*z,mu);

 dr=Dr(end);
 a.H=pack(.5150+2.883*dr+(.1695-.7018*dr)*ld,6.505-2.985*dr+(-.007969-.4299*dr)*ld,.09978+.7974*dr+(.004994-.07005*dr)*ld,.09952+.7996*dr+(.03988-.1606*dr)*ld);
 a.M=pack(44.89,.3515,.300+.4986*dr,.09981+.3710*dr+(.01998-.09041*dr)*ld);
end
names=fieldnames(a);for j=1:numel(names)
 b=a.(names{j});assert(all(b.xu>0&b.k>0&b.yu>0&b.n>=0&b.n<=1),'PISA coefficients invalid outside admissible domain.');
 assert(all(b.xu.*b.k>=b.yu*(1-1e-12)),'Conic ultimate displacement below elastic intersection.');
end
end
function a=pack(xu,k,n,yu),a=struct('xu',xu,'k',k,'n',n,'yu',yu);end
