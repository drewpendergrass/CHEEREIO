import settings_interface as si 

spc_config = si.getSpeciesConfig()

############################################################
###############CHECK RUN-IN-PLACE SETTINGS##################
############################################################

if (spc_config["DO_RUN_IN_PLACE"] == "True") and (spc_config["DO_VARON_RERUN"] == "True"):
	raise ValueError('Run-in-place and Varon et. al. rerun cannot both be turned on. Set one or both to False.')

############################################################
###############CHECK BOOLEAN CAPITALIZATION#################
############################################################

lower_case_booleans = ["DO_SPINUP",
"DO_CONTROL_RUN",
"DO_CONTROL_WITHIN_ENSEMBLE_RUNS",
"DO_ENS_SPINUP",
"ENS_SPINUP_FROM_BC_RESTART",
"SIMPLE_SCALE_FOR_FIRST_ASSIM_PERIOD",
"AMPLIFY_ENSEMBLE_SPREAD_FOR_FIRST_ASSIM_PERIOD",
"SIMPLE_SCALE_AT_END_OF_BURN_IN_PERIOD"]

for b in lower_case_booleans:
	val = spc_config[b]
	if val not in ['true','false']:
		raise ValueError(f'Setting {b} must be true or false (case sensitive); current value is {val}.')

upper_case_booleans = ["SaveLevelEdgeDiags",
"SaveStateMet",
"SaveArea",
"SaveDOFS",
"speedyCorrelationApprox",
"lognormalErrors",
"MaskCoastsGT25pctOcean",
"AV_TO_GC_GRID",
"AveragePriorAndPosterior",
"AverageScaleFactorPosteriorWithPrior",
"DO_RUN_IN_PLACE",
"DIFFERENT_RUN_IN_PLACE_FOR_BURN_IN",
"DO_VARON_RERUN",
"useLogScaleForEmissionsMaps"]

for b in upper_case_booleans:
	val = spc_config[b]
	if val not in ['True','False']:
		raise ValueError(f'Setting {b} must be True or False (case sensitive); current value is {val}.')


upper_case_boolean_dicts = ["correlatedInitialScalings",
"MaskOceanScaleFactor",
"Mask60NScaleFactor",
"Mask60SScaleFactor",
"USE_DIFFERENT_GAMMA_FOR_BURN_IN",
"Extensions"]

for b in upper_case_boolean_dicts:
	dic = spc_config[b]
	for key in dic:
		val = dic[key]
		if val not in ['True','False']:
			raise ValueError(f'Setting {b} must have True or False (case sensitive) for all values; current value for {key} within {b} is {val}.')


