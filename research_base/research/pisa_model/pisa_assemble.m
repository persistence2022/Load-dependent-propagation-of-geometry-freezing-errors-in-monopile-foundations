function a=pisa_assemble(cfg,D,t,L)
assert(D>2*t&&t>0&&L>0);ne=cfg.n_elements;assert(ne>=1&&ne==round(ne)&&cfg.gauss_order==4);
Ds=D;if ~isempty(cfg.soil_diameter_m),Ds=cfg.soil_diameter_m;assert(isscalar(Ds)&&Ds>0);end
z=linspace(0,L,ne+1)';nd=3*ne+2;nq=4*ne;
if strcmp(cfg.section,'thin_wall'),A=pi*D*t;I=pi*D^3*t/8;
elseif strcmp(cfg.section,'exact_annulus'),A=pi*(D^2-(D-2*t)^2)/4;I=pi*(D^4-(D-2*t)^4)/64;
else,error('Unknown section convention');end
EI=cfg.E_Pa*I;GA=cfg.E_Pa/(2*(1+cfg.nu))*A*cfg.shear_factor;
gp=[-.8611363115940526,-.3399810435848563,.3399810435848563,.8611363115940526];gw=[.3478548451374538,.6521451548625461,.6521451548625461,.3478548451374538];
N=sparse(nq,nd);T=N;B=N;K=sparse(nd,nd);Ks=K;w=zeros(nq,1);zq=w;
for e=1:ne
 h=z(e+1)-z(e);idx=[2*e-1,2*e,2*(ne+1)+e,2*e+1,2*e+2];
 for g=1:4
  k=(e-1)*4+g;s=(gp(g)+1)/2;[nv,nt,b]=pisa_shape(s,h);N(k,idx)=nv;T(k,idx)=nt;B(k,idx)=b;
  w(k)=gw(g)*h/2;zq(k)=z(e)+s*h;K(idx,idx)=K(idx,idx)+EI*w(k)*(b'*b);
 end
 K(idx(3),idx(3))=K(idx(3),idx(3))+GA*h;
 Ks(idx(3),idx(3))=GA*h;
end
soil=pisa_soil_at(cfg,zq);base=pisa_soil_at(cfg,L);dr=0;if cfg.model=="dunkirk",dr=soil.Dr_fraction;end
if isempty(cfg.reaction_function)
 pa=pisa_parameters(cfg.model,zq,Ds,L,dr);if cfg.model=="dunkirk",dr=base.Dr_fraction;end
 pb=pisa_parameters(cfg.model,L,Ds,L,dr);
else
 assert(~any(cfg.components(3:4)),'Custom distributed-reaction plugin requires toe reactions explicitly off.');pa=[];pb=[];
end
sc=ones(nd,1);sc(1:2:2*(ne+1))=D;S=spdiags(sc,0,nd,nd);
a=struct('z_m',z,'zq_m',zq,'weights',w,'N',N,'T',T,'B',B,'Kbeam',K,'Kshear',Ks,'EI',EI,'GA',GA,'A',A,'I',I, ...
 'soil',soil,'base_soil',base,'parameters',pa,'base_parameters',pb,'ndof',nd,'base_dofs',[2*ne+1,2*ne+2],'S',S, ...
 'soil_diameter_m',Ds,'beam_diameter_m',D,'counterfactual_soil_diameter',Ds~=D);
end
