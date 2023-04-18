import settings_interface as si 
import json

spc_config = si.getSpeciesConfig()
GC_version = int(spc_config['GC_VERSION'].split('.')[0])


############################################################
#######CHECK THAT ALL EXPECTED KEYS ARE PRESENT#############
############################################################

#Load a template that contains all expected keys
with open('../ens_config_CH4_obspack_dev.json') as f:
	template = json.load(f)

expected_keys = list(template.keys())

expected_keys = [key for key in expected_keys if 'comment' not in key] #Drop the comments, doesn't matter if they are there.

for key in expected_keys:
	if key not in spc_config:
		raise ValueError(f'Setting {key} is required in ens_config.json, even if it is not used by your simulation, but is not present.')

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
"SIMPLE_SCALE_AT_END_OF_BURN_IN_PERIOD",
"ACTIVATE_OBSPACK"]

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
"ASSIMILATE_OBS",
"AV_TO_GC_GRID",
"Extensions"]

for b in upper_case_boolean_dicts:
	dic = spc_config[b]
	for key in dic:
		val = dic[key]
		if val not in ['True','False']:
			raise ValueError(f'Setting {b} must have True or False (case sensitive) for all values; current value for key {key} within {b} is {val}.')

############################################################
###############CHECK RUN-IN-PLACE SETTINGS##################
############################################################

if (spc_config["DO_RUN_IN_PLACE"] == "True") and (spc_config["DO_VARON_RERUN"] == "True"):
	raise ValueError('Run-in-place and Varon et. al. rerun cannot both be turned on. Set one or both to False.')


############################################################
##################CHECK OBSPACK STUFF#######################
############################################################

if (GC_version==13) and (spc_config['ACTIVATE_OBSPACK'] == 'true'):
	raise ValueError('Obspack integration implemented for only GC14 and newer.')

if (not spc_config['obspack_gc_input_file'].endswith('YYYYMMDD.nc')) and (spc_config['ACTIVATE_OBSPACK'] == 'true'):
	raise ValueError(f"OBSPACK input file has to end in YYYYMMDD.nc; current noncompliant value is {spc_config['obspack_gc_input_file']}.")


