import settings_interface as si 

#USER: if you have implemented a new observation operator, add it to the operators.json file following the instructions on the Observations page in the documentation
translators = si.importObsTranslators()

obs_types = list(set(si.getSpeciesConfig()["OBS_TYPE"].values()))

for ot in obs_types:
	translator = translators[ot]()
	_ = translator.initialReadDate()