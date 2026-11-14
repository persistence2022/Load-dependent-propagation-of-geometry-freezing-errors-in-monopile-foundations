function compare_openpile_curves()
here=fileparts(mfilename('fullpath'));out=fullfile(here,'results');T=readtable(fullfile(out,'openpile_reference_curves.csv'),'TextType','string');n=height(T);ours=zeros(n,1);normerr=ours;snapped=false(n,1);
for i=1:n
 t=T(i,:);a=pisa_parameters(t.model,t.z_m,t.D_m,t.L_m,t.Dr_fraction);
 if t.model=="cowden",s=t.Su_Pa;else,s=t.sigma_Pa;end
 if t.component=="p",xs=s*t.D_m/t.G0_Pa;ys=s*t.D_m;pa=a.p;
 elseif t.component=="m",xs=s/t.G0_Pa;if t.model=="cowden",ys=s*t.D_m^2;else,ys=t.p_for_m_N_per_m*t.D_m;end;pa=a.m;
 elseif t.component=="H",xs=s*t.D_m/t.G0_Pa;ys=s*t.D_m^2;pa=a.H;
 else,xs=s/t.G0_Pa;ys=s*t.D_m^3;pa=a.M;end
 x=t.x_SI/xs;ours(i)=pisa_conic(x,pa)*ys;normerr(i)=abs(ours(i)-t.openpile_reaction_SI)/(pa.yu*ys);
 snapped(i)=abs(x-pa.xu)<.01001&&abs(x-pa.xu)>1e-5*max(pa.xu,1e-6)&&x<pa.xu;
end
T.matlab_reaction_SI=ours;T.absolute_error_over_ultimate=normerr;T.source_endpoint_snap_branch=snapped;
writetable(T,fullfile(out,'openpile_pointwise_comparison.csv'));
S=struct([]);
for model=["cowden","dunkirk"]
 for comp=["p","m","H","M"]
  j=T.model==model&T.component==comp;j2=j&~snapped;
  S=[S;struct('model',model,'component',comp,'points',sum(j),'source_snap_points',sum(j&snapped),'max_all_normalized_error',max(normerr(j)),'max_non_snap_normalized_error',max(normerr(j2)),'non_snap_passed',max(normerr(j2))<2e-5)];
 end
end
ST=struct2table(S);writetable(ST,fullfile(out,'openpile_comparison_summary.csv'));disp(ST);

assert(all(ST.non_snap_passed),'OpenPile disagreement outside documented endpoint snap requires investigation.');
end
