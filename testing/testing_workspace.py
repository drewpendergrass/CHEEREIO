import toolbox as tx

tx.makeDistTensor()

lat,lon = tx.getLatLonVals()
tx.calcDist_km(lat[0],lon[0],lat[1],lon[1])
tx.calcDist_km(lat[1],lon[1],lat[0],lon[0])

for i in range(1,9):
	for j in range(1,9):
		latlist,lonlist = tx.getLatLonList(i,j)
		for latind,lonind in zip(latlist,lonlist):
			inds = tx.getIndsOfInterest(latind,lonind)
			print(f'Inds of interest around {(latind,lonind)}, processed by core {(i,j)} are {inds}')
