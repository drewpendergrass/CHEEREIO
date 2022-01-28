import numpy as np

landmask_fraction = 0.75

gridin = {"lon":-180.0 + (0.3125*np.arange(0.0,1152.0,1.0)), "lat":-90.0 + (0.25*np.arange(0.0,721.0,1.0))}
gridout = {"lon":np.arange(-180.0,176.0, 5.0), "lat":np.concatenate([[-89.0],np.arange(-86.0,87.0, 4.0), [89.0]])}

hires_mask = np.transpose(np.genfromtxt('../templates/landmask_0p25x0p3125_gcgrid.csv',delimiter=','))

#DO 4X5
out_edge_lon = np.arange(-177.5,178.0, 5.0)
out_edge_lat = np.concatenate([[-90.0],np.arange(-88.0,89.0, 4.0)])

counts = np.zeros((len(gridout['lon']),len(gridout['lat'])))
landsum = np.zeros((len(gridout['lon']),len(gridout['lat'])))


for hireslonind,lonval in enumerate(gridin['lon']):
	if (lonval<out_edge_lon[0]) | (lonval>=out_edge_lon[-1]):
		lonind=0
	else:
		diffvals = lonval-out_edge_lon
		lonind = np.where(diffvals >= 0, diffvals, np.inf).argmin()+1
	for hireslatind,latval in enumerate(gridin['lat']):
		diffvals = latval-out_edge_lat
		latind = np.where(diffvals >= 0, diffvals, np.inf).argmin()
		counts[lonind,latind]+=1
		landsum[lonind,latind]+=hires_mask[hireslonind,hireslatind]

landfrac = landsum/counts

landfrac[landfrac>=landmask_fraction] = 1
landfrac[landfrac<landmask_fraction] = 0

landfrac = np.transpose(landfrac)
np.savetxt("../templates/landmask_4x5_gcgrid_landfraction_gt_75pct.csv",landfrac,delimiter=',')

#############################################################
#########################DO 2X2.5############################
#############################################################
out_edge_lon = np.arange(-178.75,179.0, 2.5)
out_edge_lat = np.concatenate([[-90.0],np.arange(-89.0,90.0, 2.0)])

gridout = {"lon":np.arange(-180.0,178.0, 2.5), "lat":np.concatenate([[-89.5],np.arange(-88.0,89.0, 2.0), [89.5]])}

counts = np.zeros((len(gridout['lon']),len(gridout['lat'])))
landsum = np.zeros((len(gridout['lon']),len(gridout['lat'])))

for hireslonind,lonval in enumerate(gridin['lon']):
	if (lonval<out_edge_lon[0]) | (lonval>=out_edge_lon[-1]):
		lonind=0
	else:
		diffvals = lonval-out_edge_lon
		lonind = np.where(diffvals >= 0, diffvals, np.inf).argmin()+1
	for hireslatind,latval in enumerate(gridin['lat']):
		diffvals = latval-out_edge_lat
		latind = np.where(diffvals >= 0, diffvals, np.inf).argmin()
		counts[lonind,latind]+=1
		landsum[lonind,latind]+=hires_mask[hireslonind,hireslatind]

landfrac = landsum/counts

landfrac[landfrac>=landmask_fraction] = 1
landfrac[landfrac<landmask_fraction] = 0

landfrac = np.transpose(landfrac)
np.savetxt("../templates/landmask_2x2p5_gcgrid_landfraction_gt_75pct.csv",landfrac,delimiter=',')
