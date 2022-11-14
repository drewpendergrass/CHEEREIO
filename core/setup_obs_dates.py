import settings_interface as si 

#USER: if you have implemented a new observation operator, add it to the operators.json file following the instructions on the Observations page in the documentation
translators,obsop_settings = si.importObsTranslators(return_obsop_list = True)

obs_types = list(set(si.getSpeciesConfig()["OBS_TYPE"].values()))

for ot in obs_types:
	translator = translators[ot]()
	doInitRead = obsop_settings[ot]['implements_initialReadDate'] == "True"
	if doInitRead: #Check if obs op implements initial read date function; not mandatory. 
		_ = translator.initialReadDate()