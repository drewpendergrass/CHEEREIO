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
eastchina_lon_lim = [110,135]
eastchina_lat_lim = [20,45]
southernafrica_lon_lim = [10,50]
southernafrica_lat_lim = [-35,-10]
southamerica_lon_lim = [-85,-30]
southamerica_lat_lim = [-50,15]

latlims = [conus_lat_lim,europe_lat_lim,india_lat_lim,australia_lat_lim,eastchina_lat_lim,southernafrica_lat_lim,southamerica_lat_lim]
lonlims = [conus_lon_lim,europe_lon_lim,india_lon_lim,australia_lon_lim,eastchina_lon_lim,southernafrica_lon_lim,southamerica_lon_lim]

conus_lon_ind = np.where((lon>=conus_lon_lim[0]) & (lon<=conus_lon_lim[1]))[0]
conus_lat_ind = np.where((lat>=conus_lat_lim[0]) & (lat<=conus_lat_lim[1]))[0]
europe_lon_ind = np.where((lon>=europe_lon_lim[0]) & (lon<=europe_lon_lim[1]))[0]
europe_lat_ind = np.where((lat>=europe_lat_lim[0]) & (lat<=europe_lat_lim[1]))[0]
india_lon_ind = np.where((lon>=india_lon_lim[0]) & (lon<=india_lon_lim[1]))[0]
india_lat_ind = np.where((lat>=india_lat_lim[0]) & (lat<=india_lat_lim[1]))[0]
australia_lon_ind = np.where((lon>=australia_lon_lim[0]) & (lon<=australia_lon_lim[1]))[0]
australia_lat_ind = np.where((lat>=australia_lat_lim[0]) & (lat<=australia_lat_lim[1]))[0]
eastchina_lon_ind = np.where((lon>=eastchina_lon_lim[0]) & (lon<=eastchina_lon_lim[1]))[0]
eastchina_lat_ind = np.where((lat>=eastchina_lat_lim[0]) & (lat<=eastchina_lat_lim[1]))[0]
southernafrica_lon_ind = np.where((lon>=southernafrica_lon_lim[0]) & (lon<=southernafrica_lon_lim[1]))[0]
southernafrica_lat_ind = np.where((lat>=southernafrica_lat_lim[0]) & (lat<=southernafrica_lat_lim[1]))[0]
southamerica_lon_ind = np.where((lon>=southamerica_lon_lim[0]) & (lon<=southamerica_lon_lim[1]))[0]
southamerica_lat_ind = np.where((lat>=southamerica_lat_lim[0]) & (lat<=southamerica_lat_lim[1]))[0]

latinds = [conus_lat_ind,europe_lat_ind,india_lat_ind,australia_lat_ind,eastchina_lat_ind,southernafrica_lat_ind,southamerica_lat_ind]
loninds = [conus_lon_ind,europe_lon_ind,india_lon_ind,australia_lon_ind,eastchina_lon_ind,southernafrica_lon_ind,southamerica_lon_ind]

regionnames = ['CONUS','Europe','India','Australia','EastChina','SouthernAfrica','SouthAmerica']

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

    #####Regions#########
    for latlim,lonlim,latind,lonind,rname in zip(latlims,lonlims,latinds,loninds,regionnames):
        def animate_region(i):
            daystring = timestr[i]
            titlestring = f'{variable} for {daystring}'
            plt.title(titlestring)
            temp = ensmean[i,lonind,latind]
            temp = temp[:-1, :-1] #weird old bug fix found on stackoverflow
            #mesh = m.pcolormesh(lon, lat, maptimeseries[:,:,i],latlon=True)
            mesh.set_array(temp.ravel())
            return mesh

        # call the animator.  blit=True means only re-draw the parts that have changed.
        fig = plt.figure(figsize=(10, 8))
        m = Basemap(projection='cyl', resolution='l',llcrnrlat=latlim[0], urcrnrlat=latlim[1],llcrnrlon=lonlim[0], urcrnrlon=lonlim[1])
        m.drawcountries(color='lightgray')
        m.drawcoastlines(color='lightgray')
        mesh = m.pcolormesh(lon[lonind], lat[latind], ensmean[0,lonind,latind],latlon=True,cmap=plt.cm.jet)
        plt.clim(np.min(ensmean[:,lonind,conus_lat_ind]), np.max(ensmean[:,lonind,latind]))
        plt.colorbar(label=variable);
        anim = animation.FuncAnimation(fig, animate_conus,len(time), blit=False)
        #anim = animation.FuncAnimation(fig, animate,300, blit=False) #for low memory plot
        #plt.show()

        #save as GIF

        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=anim_fps, metadata=dict(artist='Drew Pendergrass'), bitrate=800) #low res, small memory plot
        if looping:
            anim.save(f'{file_out}_{func}_{rname}.mp4', writer=writer)
        else:
            anim.save(file_out, writer=writer)



