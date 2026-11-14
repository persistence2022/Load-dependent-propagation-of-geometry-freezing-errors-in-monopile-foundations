function s=pisa_soil_at(cfg,z)
z=z(:);
if ~isempty(cfg.soil_function)
 s=cfg.soil_function(z);
else
 assert(~isempty(cfg.soil),'A documented soil profile is required.');tab=cfg.soil;
 fields={'G0_Pa','Su_Pa','sigma_v_eff_Pa','Dr_fraction'};s=struct;
 assert(min(z)>=min(tab.z_m)-1e-10&&max(z)<=max(tab.z_m)+1e-10,'Soil table does not cover embedded depth.');
 for j=1:numel(fields)
  key=fields{j};if (isstruct(tab)&&isfield(tab,key))||(istable(tab)&&ismember(key,tab.Properties.VariableNames))
   v=tab.(key);if isscalar(v),s.(key)=repmat(v,size(z));else,s.(key)=interp1(tab.z_m,v,z,'linear');end
  end
 end
end
keys=fieldnames(s);for j=1:numel(keys),v=s.(keys{j});if isscalar(v),v=v+zeros(size(z));end;s.(keys{j})=v(:);end
assert(isfield(s,'G0_Pa')&&all(isfinite(s.G0_Pa))&&all(s.G0_Pa>=0));
if cfg.model=="cowden",assert(isfield(s,'Su_Pa')&&all(s.Su_Pa>0));
else,assert(isfield(s,'sigma_v_eff_Pa')&&all(s.sigma_v_eff_Pa>=0));assert(isfield(s,'Dr_fraction'));end
end
