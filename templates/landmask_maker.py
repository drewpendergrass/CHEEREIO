import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

world=gpd.read_file('/n/home12/drewpendergrass/ne_10m_land.shp')
geom=world.geometry[0]

gridlabels = ['4.0x5.0','2.0x2.5','1x1','0.5x0.625','0.25x0.3125']
gridlabels_name = ['4x5','2x2p5','1x1','0p5x0p625','0p25x0p3125']

for i in range(2,5):
	gridlabel=gridlabels[i]
	labname=gridlabels_name[i]
	#Create traditional GEOS-Chem longitude and latitude centers, as specified by the settings in ens_config.json
	if gridlabel == '4.0x5.0':
		lon = np.arange(-180.0,176.0, 5.0)
		lat = np.concatenate([[-89.0],np.arange(-86.0,87.0, 4.0), [89.0]])
	elif gridlabel == '2.0x2.5':
		lon = np.arange(-180.0,178.0, 2.5)
		lat = np.concatenate([[-89.5],np.arange(-88.0,89.0, 2.0), [89.5]])
	elif gridlabel == '1x1':
		lon = np.arange(-179.5,180.0, 1.0)
		lat = np.arange(-89.5,90, 1.0)
	elif (gridlabel == '0.5x0.625') | (gridlabel == 'MERRA2'): #MERRA2 NATIVE GRID
		lon = -180.0 + (0.625*np.arange(0.0,576.0,1.0))
		lat = -90.0 + (0.5*np.arange(0.0,361.0,1.0))
	elif (gridlabel == 'AS_MERRA2'): #ASIA NESTED GRID FOR MERRA2
		lon = np.arange(60.0,150.01, 0.625)
		lat = np.arange(-11.0,55.01, 0.5)
	elif (gridlabel == 'EU_MERRA2'): #EUROPE NESTED GRID FOR MERRA2
		lon = np.arange(-30.0,50.01, 0.625)
		lat = np.arange(30.0,70.01, 0.5)
	elif (gridlabel == 'NA_MERRA2'): #NORTH AMERICA NESTED GRID FOR MERRA2
		lon = np.arange(-140.0,-39.99, 0.625)
		lat = np.arange(10.0,70.01, 0.5)
	elif (gridlabel == '0.25x0.3125') | (gridlabel == 'GEOSFP'):
		lon = -180.0 + (0.3125*np.arange(0.0,1152.0,1.0))
		lat = -90.0 + (0.25*np.arange(0.0,721.0,1.0))
	elif (gridlabel == 'CH_GEOSFP'): #CHINA NESTED GRID FOR GEOS-FP
		lon = np.arange(70.0,140.01, 0.3125)
		lat = np.arange(15.0,55.01, 0.25)
	elif (gridlabel == 'EU_GEOSFP'): #EU NESTED GRID FOR GEOS-FP
		lon = np.arange(-15.0,40.01, 0.3125)
		lat = np.arange(32.75,61.26, 0.25)
	elif (gridlabel == 'NA_GEOSFP'): #NA NESTED GRID FOR GEOS-FP
		lon = np.arange(-130.0,-59.99, 0.3125)
		lat = np.arange(9.75,60.01, 0.25)
	else:
		raise ValueError('Scaling factor initialization utility does not recognize grid specification.')
	lon2d, lat2d = np.meshgrid(lon,lat)
	# Reshape to 1D for easier iteration.
	lon2 = lon2d.reshape(-1)
	lat2 = lat2d.reshape(-1)
	mask = []
	# Iterate through all horizontal points in cube, and
	# check for containment within the specified geometry.
	for latval, lonval in zip(lat2, lon2):
	  this_point = Point(lonval, latval)
	  print(f"Processing {gridlabel} at ({lonval},{latval})")
	  res = geom.contains(this_point)
	  mask.append(res)
	mask = np.array(mask).reshape(lon2d.shape)
	np.savetxt(f"/n/home12/drewpendergrass/landmask_{labname}_gcgrid.csv",mask,delimiter=',')
