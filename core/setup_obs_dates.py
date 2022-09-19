import tropomi_tools as tt
import omi_tools as ot
import settings_interface as si 

obs_types = list(set(si.getSpeciesConfig()["OBS_TYPE"].values))

if "TROPOMI" in obs_types:
	translator = tt.TROPOMI_Translator()
	_ = translator.initialReadDate()

if "OMI" in obs_types:
	translator = ot.OMI_Translator()
	_ = translator.initialReadDate()
