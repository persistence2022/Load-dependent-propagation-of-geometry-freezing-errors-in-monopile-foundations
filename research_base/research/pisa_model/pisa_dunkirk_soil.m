function s=pisa_dunkirk_soil(z,Dr)

z=z(:);assert(all(z>=0));assert(isscalar(Dr));

e0=.741-(.056/.15)*(Dr-.45);sig=10090*z;pm=(1+2*.4)*sig/3;
G0=875*101300/(.3+.7*e0^2).*sqrt(pm/101300);
s=struct('G0_Pa',G0,'sigma_v_eff_Pa',sig,'Dr_fraction',Dr+zeros(size(z)));
end
