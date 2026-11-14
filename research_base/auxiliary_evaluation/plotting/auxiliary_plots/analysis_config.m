function c=analysis_config()

c.title='Load-dependent propagation of geometry-freezing errors in monopile foundations';
c.bounds=[8 11;0.06 0.14;35 60];c.names={'D_m','t_m','embedded_length_m'};
c.seeds=struct('smoke',41000,'development',41021,'calibration',41022,'test',41024,'policy',53001,'bootstrap',271828);
c.sizes=struct('smoke',6,'development',80,'calibration',60,'test',120,'policy',128);
c.loads_N_Nm=[5e6 -500e6;10e6 -300e6];
c.model_modes={'frozen_linear','updated_linear','frozen','updated'};
c.rotation_limit_rad=0.25*pi/180;c.stress_limit_Pa=345e6;c.displacement_length_ratio=1/500;
c.limit_basis='Legacy project rotation/stress thresholds; displacement=L_embedded/500 is a declared numerical-study convention, not certified site design approval.';
c.max_element_m=1;c.density_kg_m3=7850;
c.mechanism_grid=struct('D_m',linspace(8,11,13),'t_m',linspace(.06,.14,13),'L_m',[35 47.5 60]);
c.policy=struct('repetitions',8,'initial_low',12,'initial_reference',4,'budget_reference_equivalents',[8 16 24],...
 'methods',{{'Frozen','LowOnly','ReferenceOnly','Periodic','DecisionAdaptive'}},'fixed_kernel_length',0.6,'ridge',1e-6,'scale_multiplier',1.64);
c.primary_metrics={'final_reference_feasibility','conditional_feasible_mass','false_approvals_at_matched_acceptance','actual_compute_seconds'};
c.scope='New controlled static numerical data, project-manuscript soil inputs, no claim of field validation or annual failure probability.';
end
