import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap
import sys
import argparse
sys.path.append('../core')
import settings_interface as si

data = si.getSpeciesConfig()

anim_fps_emis=int(data['animation_fps_scalingfactor'])
anim_fps_conc=int(data['animation_fps_concentrations'])
pp_dir=f"{data['MY_PATH']}/{data['RUN_NAME']}/postprocess" 
concpp=f"{pp_dir}/controlvar_pp.nc"
hemcopp=f"{pp_dir}/combined_HEMCO_diagnostics.nc"
gc_version = float(data['GC_VERSION'][0:-2]) #major plus minor version
if gc_version>=14.1:
    spcconc_name = "SpeciesConcVV"
else:
    spcconc_name = "SpeciesConc" #Starting in 14.1 we have to specify VV

lognormalErrors=data['lognormalErrors']=="True"

if data['NEST'] == 'F':
    is_global = True
else:
    is_global = False

latgrid,longrid = si.getLatLonVals(data)
latgrid = np.array(latgrid)
longrid = np.array(longrid)
totallat_minmax = [np.min(latgrid),np.max(latgrid)]
totallon_minmax = [np.min(longrid),np.max(longrid)]

parser = argparse.ArgumentParser(description='Animate maps of postprocessed results; defaults to ens_config settings, but you can override.')
parser.add_argument('-func', '--statistical_function', type=str, default='all', help='Function to plot? mean, sd, max, min, range. Default: all, for all of the above.')

# create dictionary of arguments
args = parser.parse_args()
func = args.statistical_function

if is_global:
    lon = longrid
    lat = latgrid
    ## Some regional definitions.
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
    ending = ''
    length = 5
else:
    looping=False
    ending = f'_{func}.mp4'
    length = 1

variables = []
das = []
outfiles = []
times = []
timestrs = []
lons = []
lats = []
anim_fps_vals = []

#Get necessary data

#CONCENTRATIONS
ds = xr.open_dataset(concpp)
time = np.array(ds['time'])
timestr = [str(t)[0:16] for t in time]
lat = np.array(ds['lat'])
lon = np.array(ds['lon'])
anim_fps_val = anim_fps_conc
for conc_val in data['CONTROL_VECTOR_CONC']:
    variable = f'{spcconc_name}_{conc_val}'
    variables.append(variable)
    das.append(np.array(ds[variable]))
    outfiles.append(f"{pp_dir}/SpeciesConc_{conc_val}{ending}")
    times.append(time)
    timestrs.append(timestr)
    lons.append(lon)
    lats.append(lat)
    anim_fps_vals.append(anim_fps_val)

#EMISSIONS
anim_fps_val = anim_fps_emis
for emis_val in data['CONTROL_VECTOR_EMIS']:
    ds = xr.open_dataset(f"{pp_dir}/{emis_val}_SCALEFACTOR.nc")
    time = np.array(ds['time'])
    timestr = [str(t)[0:16] for t in time]
    lat = np.array(ds['lat'])
    lon = np.array(ds['lon'])
    variables.append('Scalar')
    das.append(np.array(ds['Scalar']))
    outfiles.append(f"{pp_dir}/SCALEFACTOR_{emis_val}{ending}")
    times.append(time)
    timestrs.append(timestr)
    lons.append(lon)
    lats.append(lat)
    anim_fps_vals.append(anim_fps_val)

#HEMCO
ds = xr.open_dataset(hemcopp)
time = np.array(ds['time'])
timestr = [str(t)[0:16] for t in time]
lat = np.array(ds['lat'])
lon = np.array(ds['lon'])
anim_fps_val = anim_fps_emis
for hemco_val in data['hemco_diags_to_process']:
    variables.append(hemco_val)
    das.append(np.array(ds[hemco_val]))
    outfiles.append(f"{pp_dir}/HEMCOdiag_{hemco_val}{ending}")
    times.append(time)
    timestrs.append(timestr)
    lons.append(lon)
    lats.append(lat)
    anim_fps_vals.append(anim_fps_val)


def getEnsMean(func,variable,da):
    if func=='mean':
        if lognormalErrors and (variable=='Scalar'):
            ensmean = np.exp(np.mean(np.log(da),axis=0))
        else:
            ensmean = np.mean(da,axis=0)
    elif func=='sd':
        if lognormalErrors and (variable=='Scalar'):
            ensmean = np.std(np.log(da),axis=0)
        else:
            ensmean = np.std(da,axis=0)
    elif func=='max':
        ensmean = np.max(da,axis=0)
    elif func=='min':
        ensmean = np.min(da,axis=0)
    elif func=='range':
        ensmean = np.max(da,axis=0)-np.min(da,axis=0)
    return ensmean

#####GLOBAL or full nested region########
# call the animator.  blit=True means only re-draw the parts that have changed.
m = Basemap(projection='cyl', resolution='l',llcrnrlat=totallat_minmax[0], urcrnrlat=totallat_minmax[1],llcrnrlon=totallon_minmax[0], urcrnrlon=totallon_minmax[1])

for variable,da,file_out,time,timestr,lon,lat,anim_fps in zip(variables,das,outfiles,times,timestrs,lons,lats,anim_fps_vals):

    #Loop through functions if necessary
    for i in range(length):

        fig = plt.figure(figsize=(10, 6))
        m.drawcountries(color='lightgray')
        m.drawcoastlines(color='lightgray')

        if looping:
            func = funcs[i]

        ensmean = getEnsMean(func,variable,da)

        def animate(i):
            daystring = timestr[i]
            titlestring = f'{variable} for {daystring}'
            plt.title(titlestring)
            temp = ensmean[i,:,:]
            temp = temp[:-1, :-1] #weird old bug fix found on stackoverflow
            #mesh = m.pcolormesh(lon, lat, maptimeseries[:,:,i],latlon=True)
            mesh.set_array(temp.ravel())
            return mesh

        
        #custom bwr colormap for scalings
        if variable == 'Scalar':
            cvals  = [0.0, 1.0, np.max([np.max(ensmean),1.1])]
            colors = ["blue","white","red"]
            pltnorm=plt.Normalize(min(cvals),max(cvals))
            tuples = list(zip(map(pltnorm,cvals), colors))
            cmap = LinearSegmentedColormap.from_list("", tuples)
            clim = [0.0, np.max([np.max(ensmean),1.1])]
        else:
            cmap=plt.cm.jet
            clim = [np.min(ensmean), np.max(ensmean)]

        mesh = m.pcolormesh(lon, lat, ensmean[0,:,:],latlon=True,cmap=cmap)
        plt.clim(clim[0],clim[1])
        plt.colorbar(label=variable);
        anim = animation.FuncAnimation(fig, animate,len(time), blit=False)
        #anim = animation.FuncAnimation(fig, animate,300, blit=False) #for low memory plot
        #plt.show()

        #save as GIF

        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=anim_fps, metadata=dict(artist='CHEEREIO'), bitrate=800) #low res, small memory plot
        if looping:
            anim.save(f'{file_out}_{func}.mp4', writer=writer)
        else:
            anim.save(file_out, writer=writer)

        #Close figure
        plt.close(fig)


#####Regions for global simulation#########

if is_global:
    for latlim,lonlim,latind,lonind,rname in zip(latlims,lonlims,latinds,loninds,regionnames):
        # call the animator.  blit=True means only re-draw the parts that have changed.
        m = Basemap(projection='cyl', resolution='l',llcrnrlat=latlim[0], urcrnrlat=latlim[1],llcrnrlon=lonlim[0], urcrnrlon=lonlim[1])

        for variable,da,file_out,time,timestr,lon,lat,anim_fps in zip(variables,das,outfiles,times,timestrs,lons,lats,anim_fps_vals):
            #Loop through functions if necessary
            for i in range(length):
                
                fig = plt.figure(figsize=(10, 6))
                m.drawcountries(color='lightgray')
                m.drawcoastlines(color='lightgray')

                if looping:
                    func = funcs[i]

                ensmean = getEnsMean(func,variable,da)


                def animate_region(i):
                    daystring = timestr[i]
                    titlestring = f'{variable} for {daystring}'
                    plt.title(titlestring)
                    temp = ensmean[i,latind[0]:(latind[-1]+1),lonind[0]:(lonind[-1]+1)]
                    temp = temp[:-1, :-1] #weird old bug fix found on stackoverflow
                    #mesh = m.pcolormesh(lon, lat, maptimeseries[:,:,i],latlon=True)
                    mesh.set_array(temp.ravel())
                    return mesh


                #custom bwr colormap for scalings
                if variable == 'Scalar':
                    cvals  = [0.0, 1.0, np.max([np.max(ensmean[0,latind[0]:(latind[-1]+1),lonind[0]:(lonind[-1]+1)]),1.1])]
                    colors = ["blue","white","red"]
                    pltnorm=plt.Normalize(min(cvals),max(cvals))
                    tuples = list(zip(map(pltnorm,cvals), colors))
                    cmap = LinearSegmentedColormap.from_list("", tuples)
                    clim = [0.0, np.max([np.max(ensmean[0,latind[0]:(latind[-1]+1),lonind[0]:(lonind[-1]+1)]),1.1])]
                else:
                    cmap=plt.cm.jet
                    clim = [np.min(ensmean[0,latind[0]:(latind[-1]+1),lonind[0]:(lonind[-1]+1)]), np.max(ensmean[0,latind[0]:(latind[-1]+1),lonind[0]:(lonind[-1]+1)])]
                mesh = m.pcolormesh(lon[lonind], lat[latind], ensmean[0,latind[0]:(latind[-1]+1),lonind[0]:(lonind[-1]+1)],latlon=True,cmap=cmap)
                plt.clim(clim[0],clim[1])
                plt.colorbar(label=variable);
                anim = animation.FuncAnimation(fig, animate_region,len(time), blit=False)
                #anim = animation.FuncAnimation(fig, animate,300, blit=False) #for low memory plot
                #plt.show()

                #save as GIF

                Writer = animation.writers['ffmpeg']
                writer = Writer(fps=anim_fps, metadata=dict(artist='CHEEREIO'), bitrate=800) #low res, small memory plot
                if looping:
                    anim.save(f'{file_out}_{func}_{rname}.mp4', writer=writer)
                else:
                    anim.save(file_out, writer=writer)

                #Close figure
                plt.close(fig)



