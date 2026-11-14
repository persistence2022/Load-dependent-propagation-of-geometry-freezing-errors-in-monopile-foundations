function T=field_scenarios(site)

ids=["PISA_full","PISA_p_only","PISA_G0_080","PISA_G0_120"]';
kind=repmat("PISA",4,1);gs=[1;1;.8;1.2];ss=ones(4,1);phi=nan(4,1);eps50=phi;J=phi;
if strcmpi(site,'Cowden')
 ids=[ids;"PISA_Su_080";"PISA_Su_120";"Matlock_e004_J050";"Matlock_e005_J050";"Matlock_e007_J050";"Matlock_e005_J025"];
 kind=[kind;"PISA";"PISA";repmat("Matlock",4,1)];gs=[gs;ones(6,1)];ss=[ss;.8;1.2;ones(4,1)];
 phi=[phi;nan(6,1)];eps50=[eps50;nan(2,1);.004;.005;.007;.005];J=[J;nan(2,1);.5;.5;.5;.25];
else
 ids=[ids;"API_phi32";"API_phi36";"API_phi40"];kind=[kind;repmat("API",3,1)];gs=[gs;ones(3,1)];ss=[ss;ones(3,1)];
 phi=[phi;32;36;40];eps50=[eps50;nan(3,1)];J=[J;nan(3,1)];
end
T=table(ids,kind,gs,ss,phi,eps50,J,'VariableNames',{'scenario_id','kind','G0_scale','Su_scale','phi_deg','epsilon50','J'});
T.clay_regularization=ones(height(T),1)*1e-6;
T.evidence=repmat("Predeclared parameter scenario; no response fitting",height(T),1);
T.evidence(1:2)="Published interpreted soil; explicit PISA calibration-domain extrapolation";
end
