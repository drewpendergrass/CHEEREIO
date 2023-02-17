import argparse
import numpy as np
import xarray as xr
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap
import pickle
from map_tools import *
import sys
sys.path.append('../core')
import settings_interface as si 

def str2bool(v):
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser(description='CHEEREIO runtime snapshot')
parser.add_argument("--gossipgirl", type=str2bool, nargs='?',const=True, default=False,help="Activate Gossip Girl mode.")
args = parser.parse_args()
gg_mode = args.gossipgirl

data = si.getSpeciesConfig()
pp_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/postprocess"
useControl=data['DO_CONTROL_RUN']=="true"
controlInEns = data['DO_CONTROL_WITHIN_ENSEMBLE_RUNS']=="true"
anim_fps = int(data['animation_fps_scalingfactor'])
gclat,gclon = si.getLatLonVals(data)
gclat = np.array(gclat)
gclon = np.array(gclon)

if useControl and controlInEns:
	control_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/ensemble_runs/{data['RUN_NAME']}_0000"
else:
	control_dir = None

m = Basemap(projection='cyl', resolution='l',llcrnrlat=-90, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180)

print('Generating observation and control minus observation maps...')

with open(f'{pp_dir}/SNAPSHOT_bigy_arrays_for_plotting.pkl','rb') as f:
	pickledata=pickle.load(f)

specieslist = pickledata["species"]
total_obs_in_period,total_weighted_mean_true_obs,assim_minus_obs,ctrl_minus_obs = regridBigYdata(pickledata,gclat,gclon)

print('')

for i,species in enumerate(specieslist):
	clim_abs = np.max([np.nanmax(np.abs(assim_minus_obs[i,:,:])),np.nanmax(np.abs(ctrl_minus_obs[i,:,:]))])
	plotMap(m,gclat,gclon,assim_minus_obs[i,:,:],species,f'{pp_dir}/SNAPSHOT_assim_minus_obs_{species}.png',cmap=plt.cm.seismic,clim = [-1*clim_abs,clim_abs])
	print(f'For species {species} we have, for assimilation minus observations, a mean of {np.nanmean(assim_minus_obs[i,:,:])} and a standard deviation of {np.nanstd(assim_minus_obs[i,:,:])}')
	plotMap(m,gclat,gclon,ctrl_minus_obs[i,:,:],species,f'{pp_dir}/SNAPSHOT_ctrl_minus_obs_{species}.png',cmap=plt.cm.seismic,clim = [-1*clim_abs,clim_abs]) 
	print(f'For species {species} we have, for control minus observations, a mean of {np.nanmean(ctrl_minus_obs[i,:,:])} and a standard deviation of {np.nanstd(ctrl_minus_obs[i,:,:])}')
	print('')


print('Observation and control minus observation maps complete!')
print('')

print('Starting to generate movies of scaling factors.')

path_to_sfs = glob(f'{pp_dir}/SNAPSHOT_*_SCALEFACTOR.nc')
path_to_sfs.sort()
sf_names = [pts.split('/')[-1] for pts in path_to_sfs]
sfs = {}

for name,path in zip(sf_names,path_to_sfs):
	sfs[name] = xr.open_dataset(path)



for sf in sfs:
	ds = sfs[sf]
	time = np.array(ds['time'])
	timestr = [str(t)[0:16] for t in time]
	lat = np.array(ds['lat'])
	lon = np.array(ds['lon'])
	da = np.array(ds['Scalar'])
	ensmean = np.mean(da,axis=0)

	def animate(i):
		daystring = timestr[i]
		titlestring = f'Scaling factor for {daystring}'
		plt.title(titlestring)
		temp = ensmean[i,:,:]
		temp = temp[:-1, :-1] #weird old bug fix found on stackoverflow
		#mesh = m.pcolormesh(lon, lat, maptimeseries[:,:,i],latlon=True)
		mesh.set_array(temp.ravel())
		return mesh

	#####GLOBAL#########
	# call the animator.  blit=True means only re-draw the parts that have changed.
	fig = plt.figure(figsize=(10, 6))
	m.drawcountries(color='lightgray')
	m.drawcoastlines(color='lightgray')

	#custom bwr colormap for scalings
	cvals  = [0.0, 1.0, np.max([np.max(ensmean),1.1])]
	colors = ["blue","white","red"]
	pltnorm=plt.Normalize(min(cvals),max(cvals))
	tuples = list(zip(map(pltnorm,cvals), colors))
	cmap = LinearSegmentedColormap.from_list("", tuples)
	clim = [0.0, np.max([np.max(ensmean),1.1])]

	mesh = m.pcolormesh(lon, lat, ensmean[0,:,:],latlon=True,cmap=cmap)
	plt.clim(clim[0],clim[1])
	plt.colorbar(label='Scalar');
	anim = animation.FuncAnimation(fig, animate,len(time), blit=False)

	#SAVE 
	Writer = animation.writers['ffmpeg']
	writer = Writer(fps=anim_fps, metadata=dict(artist='Drew Pendergrass'), bitrate=800) #low res, small memory plot
	sf_nodot = sf.split('.')[0]
	file_out = f'{pp_dir}/{sf_nodot}_mean_SNAPSHOT.mp4'
	anim.save(file_out, writer=writer)
	print('')
	print('')
	print(f'Saved {sf} movie of ensemble mean scalefactors out at {file_out}')
	print('')
	print('')


print('Done generating movies of scaling factors.')

if gg_mode:
	print('')
	print('')
	print('Hey Upper East Siders. Gossip Girl here.')
	print('Grab your shades and your sunblock. This one looks like a scorcher.')
	print('Our snapshot of CHEEREIO at runtime is complete. They are absolutely scandalous!')
	print(f"Don't believe me? See for yourself. View the snapshot figures and movies at {pp_dir}.")
	print('Until next time. You know you love me. XOXO —Gossip Girl.')
	print('')
else:
	print(f'Snapshot complete. Please view the figures and movies at {pp_dir}.')
