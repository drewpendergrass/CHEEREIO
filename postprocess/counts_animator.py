import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap
import pickle
import json

with open('../ens_config.json') as f:
	data = json.load(f)

pp_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/postprocess"

anim_fps = int(data['animation_fps'])

with open(f'{pp_dir}/count_arrays_for_plotting.pkl','rb') as f:
	pickledata=pickle.load(f)

arraysbase = pickledata[0:2]
dates = pickledata[2]
specieslist = pickledata[3]

with open(f"{data['MY_PATH']}/{data['RUN_NAME']}/scratch/latlon_vals.json") as f:
	ll_data = json.load(f)

gclat = np.array(ll_data['lat'])
gclon = np.array(ll_data['lon'])

filenamesbase = [f'{pp_dir}/total_raw_satellite_counts',f'{pp_dir}/total_averaged_satellite_counts']

for arrayval,filenamebase in zip(arraysbase,filenamesbase):
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
		fig = plt.figure(figsize=(10, 8))
		m = Basemap(projection='cyl', resolution='l',llcrnrlat=-90, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180)
		m.drawcountries(color='lightgray')
		m.drawcoastlines(color='lightgray')
		mesh = m.pcolormesh(gclon, gclat, arrayval[0,ind,:,:],latlon=True,cmap=plt.cm.jet)
		plt.clim(np.min(arrayval[:,ind,:,:]), np.max(arrayval[:,ind,:,:]))
		plt.colorbar(label="Count");
		anim = animation.FuncAnimation(fig, animate,len(dates), blit=False)
		Writer = animation.writers['ffmpeg']
		writer = Writer(fps=anim_fps, metadata=dict(artist='Drew Pendergrass'), bitrate=800) #low res, small memory plot
		anim.save(f'{filenamebase}_{species}.mp4', writer=writer)
