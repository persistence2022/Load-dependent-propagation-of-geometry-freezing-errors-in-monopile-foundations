function r=solve_pisa(cfg,D,t,L,loading)

tic;a=pisa_assemble(cfg,D,t,L);nd=a.ndof;ne=cfg.n_elements;
if ~isfield(loading,'control'),loading.control='force';end
dispcontrol=strcmp(loading.control,'head_displacement');
if dispcontrol
 targets=loading.w_m(:);hh=loading.eccentricity_m;unit=zeros(nd,1);unit(1)=1;unit(2)=-hh;H=zeros(size(targets));M=H;
else
 H=loading.H_N(:);M=loading.M_Nm(:);assert(numel(H)==numel(M));targets=zeros(size(H));
end
Ds=a.soil_diameter_m;
domain=struct('diameter_in_published_domain',Ds>=5&&Ds<=10,'slenderness_in_published_domain',L/Ds>=2&&L/Ds<=6,'beam_diameter_in_published_domain',D>=5&&D<=10,'density_in_published_domain',true,'load_eccentricity_in_published_domain',true);
if cfg.model=="dunkirk",domain.density_in_published_domain=all(a.soil.Dr_fraction>=.45&a.soil.Dr_fraction<=.90);end
if dispcontrol,ratio=hh/Ds;else,active=H~=0;ratio=abs(M(active)./H(active))/Ds;end
domain.load_eccentricity_in_published_domain=all(ratio>=5&ratio<=15);domain.all_in_domain=all(structfun(@(v)logical(v),domain));
domain.pisa_reaction_model=isempty(cfg.reaction_function);
domain.counterfactual_soil_diameter=a.counterfactual_soil_diameter;
if ~domain.all_in_domain&&~cfg.allow_extrapolation&&domain.pisa_reaction_model,error('PISA calibration domain exceeded. Explicit cfg.allow_extrapolation=true is required and remains flagged.');end
q=zeros(nd,1);lambda=0;Q=zeros(nd,numel(targets));iters=zeros(numel(targets),1);rels=iters;F=zeros(nd,1);free=2:nd;
sc=diag(a.S);Sf=spdiags(sc(free),0,nd-1,nd-1);
for j=1:numel(targets)
 if dispcontrol,q(1)=targets(j);else,F(1)=H(j);F(2)=M(j);end
 converged=false;
 for it=1:cfg.max_iterations
  if dispcontrol,F=lambda*unit;end
  [R,J]=pisa_residual(cfg,a,D,q,F);rr=norm(a.S'*R,inf)/max(1,norm(a.S'*F,inf));
  if rr<cfg.relative_tolerance,converged=true;break,end
  if dispcontrol
   A=[a.S'*J(:,free)*Sf,-a.S'*unit*cfg.load_scale_N];du=-(A\(a.S'*R));dq=zeros(nd,1);dq(free)=Sf*du(1:end-1);dl=du(end)*cfg.load_scale_N;
  else,dq=-a.S*((a.S'*J*a.S)\(a.S'*R));dl=0;end
  alpha=1;base_merit=norm(a.S'*R,2);
  while alpha>=cfg.line_search_min
   qt=q+alpha*dq;lt=lambda+alpha*dl;Ft=F;if dispcontrol,Ft=lt*unit;end
   Rt=pisa_residual(cfg,a,D,qt,Ft);
   if norm(a.S'*Rt,2)<=(1-1e-4*alpha)*base_merit,break,end
   alpha=alpha/2;
  end
  if alpha<cfg.line_search_min,error('Newton line search failed at path index %d, residual %.4g.',j,rr);end
  q=qt;lambda=lt;
 end
 assert(converged,'Newton failed at path index %d, residual %.4g.',j,rr);
 if dispcontrol,H(j)=lambda;M(j)=-hh*lambda;end
 Q(:,j)=q;iters(j)=it;rels(j)=rr;
end
[~,J,o]=pisa_residual(cfg,a,D,q,F);
forcebalance=sum(a.weights.*o.p_N_per_m)+o.Hb_N-H(end);
momentbalance=sum(a.weights.*(a.zq_m.*o.p_N_per_m+o.m_N))+L*o.Hb_N+o.Mb_Nm-M(end);
r=struct('response',table(H,M,Q(1,:)',Q(2,:)',iters,rels,'VariableNames',{'H_N','M_theta_Nm','head_w_m','head_theta_rad','newton_iterations','relative_residual'}), ...
 'z_m',a.z_m,'nodal_w_m',Q(1:2:2*(ne+1),:),'nodal_theta_rad',Q(2:2:2*(ne+1),:),'all_dofs',Q,'profile', ...
 table(a.zq_m,o.w_m,o.theta_rad,o.p_N_per_m,o.m_N,o.beam_moment_theta_Nm,'VariableNames',{'z_m','w_m','theta_rad','p_N_per_m','m_N','M_theta_Nm'}), ...
 'base_reactions',[o.Hb_N,o.Mb_Nm],'force_balance_error_N',forcebalance,'moment_balance_error_Nm',momentbalance, ...
 'domain',domain,'configuration',cfg,'runtime_seconds',toc,'section_A_m2',a.A,'section_I_m4',a.I,'final_tangent',J, ...
 'soil_diameter_m',Ds,'beam_diameter_m',D,'status','Executed numerical implementation; external physical validation is separate');
end
