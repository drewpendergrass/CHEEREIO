import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap

#Animate data 
def animateData(m,data,file_out,lon,lat,anim_fps = 8, variable = 'Variable',timestr = None, cmin = None, cmax=None, bwr_cmap=False):

    fig = plt.figure(figsize=(10, 6))
    m.drawcountries(color='lightgray')
    m.drawcoastlines(color='lightgray')

    def animate(i):
    	if timestr is not None:
        	daystring = timestr[i]
        	titlestring = f'{variable} for {daystring}'
        else:
        	titlestring = variable
        plt.title(titlestring)
        temp = data[i,:,:]
        temp = temp[:-1, :-1] #weird old bug fix found on stackoverflow
        #mesh = m.pcolormesh(lon, lat, maptimeseries[:,:,i],latlon=True)
        mesh.set_array(temp.ravel())
        return mesh

    #custom bwr colormap for scalings
    if bwr_cmap:
    	if cmax is None:
    		cmax = np.nanmax([np.nanmax(data),1.1])
        cvals  = [0.0, 1.0, cmax]
        colors = ["blue","white","red"]
        pltnorm=plt.Normalize(min(cvals),max(cvals))
        tuples = list(zip(map(pltnorm,cvals), colors))
        cmap = LinearSegmentedColormap.from_list("", tuples)
        clim = [0.0, cmax]
    else:
        cmap=plt.cm.jet
        if cmin is None:
        	cmin = np.nanmin(data)
        if cmax is None:
        	cmax = np.nanmax(data)
        clim = [cmin, cmax]

    mesh = m.pcolormesh(lon, lat, data[0,:,:],latlon=True,cmap=cmap)
    plt.clim(clim[0],clim[1])
    plt.colorbar(label=variable);
    anim = animation.FuncAnimation(fig, animate,np.shape(data)[0], blit=False)
    #anim = animation.FuncAnimation(fig, animate,300, blit=False) #for low memory plot
    #plt.show()

    #save as GIF

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=anim_fps, metadata=dict(artist='CHEEREIO'), bitrate=800) #low res, small memory plot
    anim.save(file_out, writer=writer)

    #Close figure
    plt.close(fig)

