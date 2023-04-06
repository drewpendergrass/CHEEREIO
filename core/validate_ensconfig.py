import settings_interface as si 

spc_config = si.getSpeciesConfig()

if (spc_config["DO_RUN_IN_PLACE"] == "True") and (spc_config["DO_VARON_RERUN"] == "True"):
	raise ValueError('Run-in-place and Varon et. al. rerun cannot both be turned on. Set one or both to False.')
	
