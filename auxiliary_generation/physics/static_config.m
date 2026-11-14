function cfg=static_config()

cfg=struct('E_Pa',210e9,'rho_kg_m3',7850,'D0_m',9,'H_N',5e6,'M_Nm',-500e6,...
 'max_element_m',1,'clay_regularization',1e-4,'load_steps',12,'newton_tol',1e-8,...
 'max_iterations',80,'soil_model','layered_py','linear_k_N_m2',54e6,...
 'tip_boundary','free','return_profile',true,'py_linear',false,'linearization_fraction',0.5);
cfg.layer_bottom_m=cumsum([1.5 6.7 5.6 8.6 7.5 11.2 28.9]);
cfg.layer_type=[1 2 2 2 1 1 2];
cfg.gamma_N_m3=[7.0 8.4 9.1 9.2 9.2 9.5 9.7]*1e3;
cfg.cu_Pa=[10 0 0 0 60 84 0]*1e3;
cfg.e50=[0.022 0 0 0 0.010 0.010 0];
cfg.phi_deg=[0 28 30 31 0 0 33];
cfg.k_N_m3=[0 7800 8100 8400 0 0 9000]*1e3;
cfg.provenance='Project Chinese manuscript Table 2: explicitly reported units; manuscript-defined study profile, underlying geotechnical report not available. Not FEA/input.json.';
cfg.scope='Euler-Bernoulli elastic beam, monotonic lateral p-y skeleton, no axial/lateral coupling, dynamic feedback, hysteresis or fatigue';
end
