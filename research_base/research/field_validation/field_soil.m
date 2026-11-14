function s=field_soil(z,site,P,scenario)

z=z(:);s=struct;s.G0_Pa=interp_profile(P.Gz,P.G,z)*1e6*scenario.G0_scale;
if strcmpi(site,'Dunkirk')
 s.Su_Pa=zeros(size(z));s.Dr_fraction=.75+.25*(z<3);
 s.sigma_v_eff_Pa=17.1e3*min(z,5.4)+(19.9-9.81)*1e3*max(z-5.4,0);
else
 s.Su_Pa=interp_profile(P.Sz,P.S,z)*1e3*scenario.Su_scale;s.Dr_fraction=zeros(size(z));
 s.sigma_v_eff_Pa=21.19e3*min(z,1)+(21.19-9.81)*1e3*max(z-1,0);
end
end
function v=interp_profile(zs,vs,z)
assert(max(z)<=max(zs)+1e-10,'No deep extrapolation of published soil profile.');

cuts=[0;find(diff(zs)==0);numel(zs)];v=nan(size(z));zz=max(z,min(zs));
for j=1:numel(cuts)-1
 ix=cuts(j)+1:cuts(j+1);mask=zz>=zs(ix(1))&zz<=zs(ix(end));
 v(mask)=interp1(zs(ix),vs(ix),zz(mask),'linear');
end
assert(all(isfinite(v)&v>0));
end
