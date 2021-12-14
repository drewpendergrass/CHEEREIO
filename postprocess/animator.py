import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap
import sys

file_in=str(sys.argv[1])
#Scalar for scalefactor files, SpeciesConc_SPEC for con files
variable=str(sys.argv[2])
func=str(sys.argv[3])
#if doing all, leave off ending; will append _FUNCTION.mp4
file_out=str(sys.argv[4])
anim_fps=int(sys.argv[5])
ds = xr.open_dataset(file_in)
time = np.array(ds['time'])
timestr = [str(t)[0:16] for t in time]
lat = np.array(ds['lat'])
lon = np.array(ds['lon'])
da = np.array(ds[variable])

if func == 'all':
    looping = True
    funcs = ['mean','sd','max','min','range']
    length = 5
else:
    looping=False
    length = 1

for i in range(length):
    if looping:
        func = funcs[i]
    if func=='mean':
    	ensmean = np.mean(da,axis=0)
    elif func=='sd':
    	ensmean = np.std(da,axis=0)
    elif func=='max':
    	ensmean = np.max(da,axis=0)
    elif func=='min':
    	ensmean = np.min(da,axis=0)
    elif func=='range':
        ensmean = np.max(da,axis=0)-np.min(da,axis=0)

    def animate(i):
        daystring = timestr[i]
        titlestring = f'{variable} for {daystring}'
        plt.title(titlestring)
        temp = ensmean[i,:,:]
        temp = temp[:-1, :-1] #weird old bug fix found on stackoverflow
        #mesh = m.pcolormesh(lon, lat, maptimeseries[:,:,i],latlon=True)
        mesh.set_array(temp.ravel())
        return mesh

    # call the animator.  blit=True means only re-draw the parts that have changed.
    fig = plt.figure(figsize=(10, 8))
    m = Basemap(projection='cyl', resolution='l',llcrnrlat=-90, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180)
    m.drawcountries(color='lightgray')
    m.drawcoastlines(color='lightgray')
    mesh = m.pcolormesh(lon, lat, ensmean[0,:,:],latlon=True,cmap=plt.cm.jet)
    plt.clim(np.min(ensmean), np.max(ensmean))
    plt.colorbar(label=variable);
    anim = animation.FuncAnimation(fig, animate,len(time), blit=False)
    #anim = animation.FuncAnimation(fig, animate,300, blit=False) #for low memory plot
    #plt.show()

    #save as GIF

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=anim_fps, metadata=dict(artist='Drew Pendergrass'), bitrate=800) #low res, small memory plot
    if looping:
        anim.save(f'{file_out}_{func}.mp4', writer=writer)
    else:
        anim.save(file_out, writer=writer)


