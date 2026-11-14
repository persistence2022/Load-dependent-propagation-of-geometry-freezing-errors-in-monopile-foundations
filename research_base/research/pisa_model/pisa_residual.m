function [R,J,o]=pisa_residual(cfg,a,D,q,F)
D=a.soil_diameter_m;
if isempty(cfg.reaction_function)
 [p,m,dp,dmt,dmw]=pisa_reactions(cfg.model,a.N*q,a.T*q,D,a.soil,a.parameters);
else
 [p,m,dp,dmt,dmw]=cfg.reaction_function(a.zq_m,a.N*q,a.T*q,D,a.soil);
 assert(all(isfinite([p(:);m(:);dp(:);dmt(:);dmw(:)]))&&isequal(size(p),size(a.zq_m)),'Custom reaction dimensions/values invalid.');
end
cm=cfg.components;p=p*cm(1);dp=dp*cm(1);m=m*cm(2);dmt=dmt*cm(2);dmw=dmw*cm(2);
w=a.weights;nq=numel(w);

curvature=pisa_beam_curvature(a,q);
R=a.B'*(w.*(a.EI*curvature))+a.Kshear*q+a.N'*(w.*p)+a.T'*(w.*m)-F;
J=a.Kbeam+a.N'*spdiags(w.*dp,0,nq,nq)*a.N+a.T'*spdiags(w.*dmt,0,nq,nq)*a.T+a.T'*spdiags(w.*dmw,0,nq,nq)*a.N;
id=a.base_dofs;Hb=0;Mb=0;kh=0;km=0;
if any(cm(3:4))
 s=a.base_soil;if cfg.model=="cowden",stress=s.Su_Pa;else,stress=s.sigma_v_eff_Pa;end
 vb=stress*D/s.G0_Pa;tb=stress/s.G0_Pa;
 [hb,kh]=pisa_conic(q(id(1))/vb,a.base_parameters.H);[mb,km]=pisa_conic(q(id(2))/tb,a.base_parameters.M);
 Hb=hb*stress*D^2*cm(3);Mb=mb*stress*D^3*cm(4);kh=kh*s.G0_Pa*D*cm(3);km=km*s.G0_Pa*D^3*cm(4);
end
R(id)=R(id)+[Hb;Mb];J(id,id)=J(id,id)+diag([kh,km]);
if nargout>2,o=struct('p_N_per_m',p,'m_N',m,'Hb_N',Hb,'Mb_Nm',Mb,'w_m',a.N*q,'theta_rad',a.T*q,'beam_moment_theta_Nm',a.EI*curvature);end
end
