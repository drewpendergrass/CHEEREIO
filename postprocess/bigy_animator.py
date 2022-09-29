import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap
import pickle
import sys 
sys.path.append('../core')
import settings_interface as si 

data = si.getSpeciesConfig()
gclat,gclon = si.getLatLonVals(data)
gclat = np.array(gclat)
gclon = np.array(gclon)
pp_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/postprocess"

anim_fps = int(data['animation_fps_scalingfactor'])

#Saving albedo isn't a default option for all runs, so have to check if it is even in the ens config.
if "postprocess_save_albedo" in data:
	postprocess_save_albedo = data['postprocess_save_albedo']=="True"
else:
	postprocess_save_albedo = False

with open(f'{pp_dir}/bigy_arrays_for_plotting.pkl','rb') as f:
	pickledata=pickle.load(f)

dates = pickledata["dates"]
specieslist = pickledata["species"]
total_obs_count=pickledata["obscount"]
total_averaged_obs=pickledata["obscount_avg"]
true_obs = pickledata["obs"]
sim_obs = pickledata["sim_obs"]

arraysbase=[total_obs_count,total_averaged_obs,true_obs,sim_obs]
filenamesbase = [f'{pp_dir}/total_raw_observation_counts',f'{pp_dir}/total_aggregated_observation_counts',f'{pp_dir}/actual_observations',f'{pp_dir}/simulated_observations']
labelnames = ['Count','Count','OBSUNIT', 'OBSUNIT']

if postprocess_save_albedo:
	total_swir = pickledata["swir_albedo"]
	arraysbase.append(total_swir)
	filenamesbase.append('{pp_dir}/averaged_albedo_SWIR')
	labelnames.append('Albedo')
	total_nir = pickledata["nir_albedo"]
	arraysbase.append(total_nir)
	filenamesbase.append('{pp_dir}/averaged_albedo_NIR')
	labelnames.append('Albedo')
	total_blended = pickledata["blended_albedo"]
	arraysbase.append(total_blended)
	filenamesbase.append('{pp_dir}/averaged_blended_albedo')
	labelnames.append('Albedo')


for arrayval,filenamebase,labelname in zip(arraysbase,filenamesbase,labelnames):
	for ind,species in enumerate(specieslist):
		def animate(i):
			daystring = dates[i]
			titlestring = f'{daystring}'
			plt.title(titlestring)
			temp = arrayval[i,ind,:,:]
			temp = temp[:-1, :-1] #weird old bug fix found on stackoverflow
			#mesh = m.pcolormesh(lon, lat, maptimeseries[:,:,i],latlon=True)
			mesh.set_array(temp.ravel())
			return mesh
		# call the animator.  blit=True means only re-draw the parts that have changed.
		fig = plt.figure(figsize=(10, 6))
		m = Basemap(projection='cyl', resolution='l',llcrnrlat=-90, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180)
		m.drawcountries(color='lightgray')
		m.drawcoastlines(color='lightgray')
		mesh = m.pcolormesh(gclon, gclat, arrayval[0,ind,:,:],latlon=True,cmap=plt.cm.jet)
		plt.clim(np.min(arrayval[:,ind,:,:]), np.max(arrayval[:,ind,:,:]))
		if labelname == 'OBSUNIT':
			label = f"{data['OBSERVED_SPECIES'][species]} ({data['OBSERVATION_UNITS'][species]})"
		else:
			label = labelname
		plt.colorbar(label=label);
		anim = animation.FuncAnimation(fig, animate,len(dates), blit=False)
		Writer = animation.writers['ffmpeg']
		writer = Writer(fps=anim_fps, metadata=dict(artist='CHEEREIO'), bitrate=800) #low res, small memory plot
		anim.save(f'{filenamebase}_{species}.mp4', writer=writer)
