import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm,LinearSegmentedColormap
import pickle
from glob import glob

def plotMap(m,lat,lon,flat,labelname,outfile,clim=None,cmap=None,useLog=False):
	fig = plt.figure(figsize=(10, 6))
	m.drawcountries(color='lightgray')
	m.drawcoastlines(color='lightgray')
	if cmap is None:
		cmap = plt.cm.jet
	if useLog:
		flat[flat<=0] = np.nan
		mesh = m.pcolormesh(lon, lat, flat,latlon=True,cmap=cmap,norm=LogNorm())
	else:
		mesh = m.pcolormesh(lon, lat, flat,latlon=True,cmap=cmap)
	if clim is not None:
		plt.clim(clim[0],clim[1])
	else:
		plt.clim(np.nanmin(flat), np.nanmax(flat))
	plt.colorbar(label=labelname)
	fig.savefig(outfile)

def plotScaleFactor(m,lat,lon,ppdir, plotMonthStartOnly=True):
	files = glob(f'{ppdir}/*_SCALEFACTOR.nc')
	files.sort()
	sf_names = [pts.split('/')[-1][0:-15] for pts in files]
	for file,name in zip(files,sf_names):
		ds = xr.load_dataset(file)
		dates = ds['time'].values
		scalar = ds['Scalar'].values
		if plotMonthStartOnly:
			years = dates.astype('datetime64[Y]').astype(int) + 1970
			months = dates.astype('datetime64[M]').astype(int) % 12 + 1
			ym = (years*100) + months
			_, ind = np.unique(ym,return_index = True)
			dates = dates[ind]
			scalar = scalar[:,ind,:,:]
		scalar = np.mean(scalar,axis=0) #average across ensemble
		timelabels = [str(timeval)[0:13] for timeval in dates]
		#Make custom blue-white-red colorbar centered at one
		cvals  = [0.0, 1.0, np.max([np.max(scalar),1.1])]
		colors = ["blue","white","red"]
		pltnorm=plt.Normalize(min(cvals),max(cvals))
		tuples = list(zip(map(pltnorm,cvals), colors))
		cmap = LinearSegmentedColormap.from_list("", tuples)
		for i,dateval in enumerate(timelabels):
			plotMap(m,lat,lon,scalar[i,:,:],'Scaling factor',f'{ppdir}/{name}_{dateval}_scalefactor.png',clim=[0,np.max([np.max(scalar),1.1])],cmap=cmap)

