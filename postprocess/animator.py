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

conus_lon_lim = [-130,-65]
conus_lat_lim = [20,50]
europe_lon_lim = [-10,40]
europe_lat_lim = [40,65]
india_lon_lim = [65,95]
india_lat_lim = [5,30]
australia_lon_lim = [110,155]
australia_lat_lim = [-40,-10]

conus_lon_ind = np.where((lon>=conus_lon_lim[0]) & (lon<=conus_lon_lim[1]))[0]
conus_lat_ind = np.where((lat>=conus_lat_lim[0]) & (lat<=conus_lat_lim[1]))[0]

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

    #####GLOBAL#########
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

    #####CONUS+#########
    def animate_conus(i):
        daystring = timestr[i]
        titlestring = f'{variable} for {daystring}'
        plt.title(titlestring)
        temp = ensmean[i,conus_lon_ind,conus_lat_ind]
        temp = temp[:-1, :-1] #weird old bug fix found on stackoverflow
        #mesh = m.pcolormesh(lon, lat, maptimeseries[:,:,i],latlon=True)
        mesh.set_array(temp.ravel())
        return mesh

    # call the animator.  blit=True means only re-draw the parts that have changed.
    fig = plt.figure(figsize=(10, 8))
    m = Basemap(projection='cyl', resolution='l',llcrnrlat=20, urcrnrlat=50,llcrnrlon=-130, urcrnrlon=-65)
    m.drawcountries(color='lightgray')
    m.drawcoastlines(color='lightgray')
    mesh = m.pcolormesh(lon[conus_lon_ind], lat[conus_lat_ind], ensmean[0,conus_lon_ind,conus_lat_ind],latlon=True,cmap=plt.cm.jet)
    plt.clim(np.min(ensmean[:,conus_lon_ind,conus_lat_ind]), np.max(ensmean[:,conus_lon_ind,conus_lat_ind]))
    plt.colorbar(label=variable);
    anim = animation.FuncAnimation(fig, animate_conus,len(time), blit=False)
    #anim = animation.FuncAnimation(fig, animate,300, blit=False) #for low memory plot
    #plt.show()

    #save as GIF

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=anim_fps, metadata=dict(artist='Drew Pendergrass'), bitrate=800) #low res, small memory plot
    if looping:
        anim.save(f'{file_out}_{func}_CONUS.mp4', writer=writer)
    else:
        anim.save(file_out, writer=writer)



