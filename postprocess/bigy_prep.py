import pandas as pd
import numpy as np
import pickle
import postprocess_tools as pt
import sys 
sys.path.append('../core')
import settings_interface as si 

data = si.getSpeciesConfig()

pp_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/postprocess"

with open(f"{pp_dir}/bigY.pkl",'rb') as f:
	bigy=pickle.load(f)

gclat,gclon = si.getLatLonVals(data)
gclat = np.array(gclat)
gclon = np.array(gclon)

#Saving albedo isn't a default option for all runs, so have to check if it is even in the ens config.
if "postprocess_save_albedo" in data:
	postprocess_save_albedo = data['postprocess_save_albedo']=="True"
else:
	postprocess_save_albedo = False

useControl=data['DO_CONTROL_RUN']=="true"
nEnsemble = int(data['nEnsemble'])

arraysbase = pt.makeBigYArrays(bigy,gclat,gclon,nEnsemble,postprocess_save_albedo=postprocess_save_albedo,useControl=useControl)

f = open(f'{pp_dir}/bigy_arrays_for_plotting.pkl',"wb")
pickle.dump(arraysbase,f)
f.close()
