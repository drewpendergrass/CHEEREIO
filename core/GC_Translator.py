import numpy as np
import xarray as xr
from glob import glob
import toolbox as tx 
import settings_interface as si 
from datetime import date,datetime,timedelta

#This class contains useful methods for getting data from GEOS-Chem restart files and 
#emissions scaling factor netCDFs. After initialization it contains the necessary data
#and can output it in useful ways to other functions in the LETKF procedure.
class GC_Translator(object):
	def __init__(self, path_to_rundir,timestamp,computeStateVec = False,verbose=1):
		#self.latinds,self.loninds = si.getLatLonList(ensnum)
		self.filename = f'{path_to_rundir}GEOSChem.Restart.{timestamp}z.nc4'
		self.timestamp=timestamp
		self.timestamp_as_date = np.datetime64(f'{timestamp[0:4]}-{timestamp[4:6]}-{timestamp[6:8]}T{timestamp[9:11]}:{timestamp[11:13]}:00')
		self.timestring = f'minutes since {timestamp[0:4]}-{timestamp[4:6]}-{timestamp[6:8]} {timestamp[9:11]}:{timestamp[11:13]}:00'
		self.restart_ds = xr.load_dataset(self.filename)
		self.emis_sf_filenames = glob(f'{path_to_rundir}*_SCALEFACTOR.nc')
		self.verbose=verbose
		if self.verbose>=3:
			self.num = path_to_rundir.split('_')[-1][0:4]
			print(f"GC_translator number {self.num} has been called for directory {path_to_rundir} and restart {self.filename}; construction beginning")
		self.emis_ds_list = {}
		for file in self.emis_sf_filenames:
			name = '_'.join(file.split('/')[-1].split('_')[0:-1])
			self.emis_ds_list[name] = xr.load_dataset(file)
			if self.verbose>=3:
				print(f"GC_translator number {self.num} has loaded scaling factors for {name}")
		if computeStateVec:
			self.buildStateVector()
		else:
			self.statevec = None
			self.statevec_lengths = None #Until state vector is initialized this variable is None
		if self.verbose>=3:
			print(f"GC_Translator number {self.num} construction complete.")
	#Since only one timestamp, returns in format lev,lat,lon
	def getSpecies3Dconc(self, species):
		da = np.array(self.restart_ds[f'SpeciesRst_{species}']).squeeze()
		if self.verbose>=3:
			print(f"GC_Translator number {self.num} got 3D conc for species {species} which are of dimension {np.shape(da)}.")
		return da
	def setSpecies3Dconc(self, species, conc3d):
		baseshape = np.shape(conc3d)
		conc4d = conc3d.reshape(np.concatenate([np.array([1]),baseshape]))
		if self.verbose>=3:
			print(f"GC_Translator number {self.num} set 3D conc for species {species} which are of dimension {np.shape(conc4d)}.")
		self.restart_ds[f'SpeciesRst_{species}'] = (["time","lev","lat","lon"],conc4d,{"long_name":f"Dry mixing ratio of species {species}","units":"mol mol-1 dry","averaging_method":"instantaneous"})
	def getLat(self):
		return np.array(self.restart_ds['lat'])
	def getLon(self):
		return np.array(self.restart_ds['lon'])
	def getLev(self):
		return np.array(self.restart_ds['lev'])
	def getRestartTime(self):
		return np.array(self.restart_ds['time'])
	def getEmisTime(self):
		return np.array(list(self.emis_ds_list.values())[0]['time'])
	#Get the emissions from the timestamp nearest to the one supplied by the user.
	def getEmisSF(self, species):
		da = self.emis_ds_list[species]['Scalar']
		time_array = self.getEmisTime()
		ind_closest = np.argmin(np.abs(time_array-self.timestamp_as_date))
		return np.array(da)[ind_closest,:,:].squeeze()
	def getEmisLat(self, species):
		return np.array(self.emis_ds_list[species]['lat'])
	def getEmisLon(self, species):
		return np.array(self.emis_ds_list[species]['lon'])
	#Add 2d emissions scaling factors to the end of the emissions scaling factor
	def addEmisSF(self, species, emis2d, assim_time):
		timelist = self.getEmisTime()
		last_time = timelist[-1]
		#new_last_time = last_time+np.timedelta64(assim_time,'h') #Add assim time hours to the last timestamp
		tstr = f'{self.timestamp[0:4]}-{self.timestamp[4:6]}-{self.timestamp[6:8]}T{self.timestamp[9:11]}:{self.timestamp[11:13]}:00.000000000'
		new_last_time = np.datetime64(tstr)
		if si.getSpeciesConfig()['DO_ENS_SPINUP']=='true':
			START_DATE = si.getSpeciesConfig()['ENS_SPINUP_START']
		else:
			START_DATE = si.getSpeciesConfig()['START_DATE']
		orig_timestamp = f'{START_DATE[0:4]}-{START_DATE[4:6]}-{START_DATE[6:8]}' #Start date from  JSON
		END_DATE = si.getSpeciesConfig()['END_DATE']
		end_timestamp = f'{END_DATE[0:4]}-{END_DATE[4:6]}-{END_DATE[6:8]}'
		#Create dataset with this timestep's scaling factors
		ds = xr.Dataset(
			{"Scalar": (("time","lat","lon"), np.expand_dims(emis2d,axis = 0),{"long_name": "Scaling factor", "units":"1"})},
			coords={
				"time": (["time"], np.array([new_last_time]), {"long_name": "time", "calendar": "standard", "units":f"hours since {orig_timestamp} 00:00:00"}),
				"lat": (["lat"], self.getEmisLat(species),{"long_name": "Latitude", "units":"degrees_north"}),
				"lon": (["lon"], self.getEmisLon(species),{"long_name": "Longitude", "units":"degrees_east"})
			},
			attrs={
				"Title":"CHEEREIO scaling factors",
				"Conventions":"COARDS",
				"Format":"NetCDF-4",
				"Model":"GENERIC",
				"NLayers":"1",
				"History":f"The LETKF utility added new scaling factors on {str(date.today())}",
				"Start_Date":f"{orig_timestamp}",
				"Start_Time":"0",
				"End_Date":f"{end_timestamp}",
				"End_Time":"0"
			}
		)
		self.emis_ds_list[species] = xr.concat([self.emis_ds_list[species],ds],dim = 'time') #Concatenate
	def buildStateVector(self):
		if self.verbose>=3:
			print("*****************************************************************")
			print(f"GC_Translator number {self.num} is starting build of statevector!")
		species_config = si.getSpeciesConfig()
		statevec_components = []
		for spec_conc in species_config['STATE_VECTOR_CONC']:
			statevec_components.append(self.getSpecies3Dconc(spec_conc).flatten())
		#If no scaling factor files, append 1s because this is a nature directory
		if len(self.emis_sf_filenames)==0:
			lenones = len(self.getLat())*len(self.getLon())*len(species_config['CONTROL_VECTOR_EMIS'])
			statevec_components.append(np.ones(lenones))
		else:
			for spec_emis in species_config['CONTROL_VECTOR_EMIS'].keys():
				statevec_components.append(self.getEmisSF(spec_emis).flatten())
		self.statevec_lengths = np.array([len(vec) for vec in statevec_components])
		self.statevec = np.concatenate(statevec_components)
		if self.verbose>=3:
			print(f"GC_Translator number {self.num} has built statevector; it is of dimension {np.shape(self.statevec)}.")
			print("*****************************************************************")
	def getLocalizedStateVectorIndices(self,latind,lonind):
		surr_latinds, surr_loninds = tx.getIndsOfInterest(latind,lonind)
		if self.verbose>=3:
			print(f"GC_Translator is getting localized statevec indices surrounding {(latind,lonind)} (lat/lon inds have shapes {np.shape(surr_latinds)}/{np.shape(surr_loninds)}); Lat inds are {surr_latinds} and lon inds are {surr_loninds}.")
		levcount = len(self.getLev())
		latcount = len(self.getLat())
		loncount = len(self.getLon())
		totalcount = levcount*latcount*loncount
		dummy3d = np.arange(0, totalcount).reshape((levcount,latcount,loncount))
		dummywhere_flat = dummy3d[:,surr_latinds,surr_loninds].flatten()
		if self.verbose>=3:
			print(f"Within a flattened 3D dummy cube, {len(dummywhere_flat)} entries are valid.")
		dummy2d = np.arange(0, latcount*loncount).reshape((latcount,loncount))
		dummy2dwhere_flat = dummy2d[surr_latinds,surr_loninds].flatten()
		if self.verbose>=3:
			print(f"Within a flattened 2D dummy square, {len(dummy2dwhere_flat)} entries are valid.")
		species_config = si.getSpeciesConfig()
		conccount = len(species_config['STATE_VECTOR_CONC'])
		emcount = len(species_config['CONTROL_VECTOR_EMIS'])
		ind_collector = []
		cur_offset = 0
		for i in range(conccount):
			ind_collector.append((dummywhere_flat+cur_offset))
			cur_offset+=totalcount
		for i in range(emcount):
			ind_collector.append((dummy2dwhere_flat+cur_offset))
			cur_offset+=(latcount*loncount)
		statevecinds = np.concatenate(ind_collector)
		if self.verbose>=3:
			print(f"There are a total of {len(statevecinds)}/{len(self.statevec)} selected from total statevec.")
		return statevecinds
	def getColumnIndicesFromFullStateVector(self,latind,lonind):
		if self.verbose>=3:
			print(f"GC_Translator is getting column statevec indices FOR FULL VECTOR at {(latind,lonind)}.")
		levcount = len(self.getLev())
		latcount = len(self.getLat())
		loncount = len(self.getLon())
		totalcount = levcount*latcount*loncount
		dummy3d = np.arange(0, totalcount).reshape((levcount,latcount,loncount))
		dummywhere_flat = dummy3d[:,latind,lonind].flatten()
		if self.verbose>=3:
			print(f"Within a flattened 3D dummy cube, {len(dummywhere_flat)} entries are valid.")
		dummy2d = np.arange(0, latcount*loncount).reshape((latcount,loncount))
		dummy2dwhere_flat = dummy2d[latind,lonind]
		if self.verbose>=3:
			print(f"Within a flattened 2D dummy square, {dummy2dwhere_flat} is sole valid entry.")
		species_config = si.getSpeciesConfig()
		conccount = len(species_config['STATE_VECTOR_CONC'])
		emcount = len(species_config['CONTROL_VECTOR_EMIS'])
		ind_collector = []
		cur_offset = 0
		for i in range(conccount):
			ind_collector.append(dummywhere_flat+cur_offset)
			cur_offset+=totalcount
		for i in range(emcount):
			ind_collector.append(np.array([dummy2dwhere_flat+cur_offset]))
			cur_offset+=(latcount*loncount)
		statevecinds = np.concatenate(ind_collector)
		if self.verbose>=3:
			print(f"There are a total of {len(statevecinds)}/{len(self.statevec)} selected from total statevec.")
		return statevecinds
	def getSpeciesConcIndicesInColumn(self,species):
		levcount = len(self.getLev())
		species_config = si.getSpeciesConfig()
		cur_offset = 0
		for ind,spec in enumerate(species_config['STATE_VECTOR_CONC']):
			if species == spec:
				return np.arange(cur_offset,cur_offset+levcount)
			cur_offset+=levcount
		return None #If loop doesn't terminate we did not find the species
	def getSpeciesEmisIndicesInColumn(self,species):
		levcount = len(self.getLev())
		species_config = si.getSpeciesConfig()
		cur_offset = len(species_config['STATE_VECTOR_CONC'])*levcount
		for ind,spec in enumerate(species_config['CONTROL_VECTOR_EMIS']):
			if species == spec:
				return cur_offset
			cur_offset+=1
		return None #If loop doesn't terminate we did not find the species
	def getColumnIndicesFromLocalizedStateVector(self,latind,lonind):
		surr_latinds, surr_loninds = tx.getIndsOfInterest(latind,lonind)
		if self.verbose>=3:
			print(f"GC_Translator is getting column statevec indices surrounding {(latind,lonind)} (lat/lon inds have shapes {np.shape(surr_latinds)}/{np.shape(surr_loninds)}); Lat inds are {surr_latinds} and lon inds are {surr_loninds}.")
		levcount = len(self.getLev())
		latcount = len(self.getLat())
		loncount = len(self.getLon())
		totalcount = levcount*latcount*loncount
		dummy3d = np.arange(0, totalcount).reshape((levcount,latcount,loncount))
		dummywhere_flat = dummy3d[:,surr_latinds,surr_loninds].flatten()
		dummywhere_flat_column = dummy3d[:,latind,lonind].flatten()
		dummywhere_match = np.where(np.in1d(dummywhere_flat,dummywhere_flat_column))[0]
		if self.verbose>=3:
			print(f"Within a flattened 3D dummy cube, {len(dummywhere_flat_column)} entries are valid in the column.")
			print(f"Matched {len(dummywhere_match)} entries in the overall flattened and subsetted column; values are {dummywhere_match}")
		dummy2d = np.arange(0, latcount*loncount).reshape((latcount,loncount))
		dummy2dwhere_flat = dummy2d[surr_latinds,surr_loninds].flatten()
		dummy2dwhere_flat_column = dummy2d[latind,lonind]
		dummy2dwhere_match = np.where(np.in1d(dummy2dwhere_flat,dummy2dwhere_flat_column))[0]
		if self.verbose>=3:
			print(f"Within a flattened 2D dummy square, {dummy2dwhere_flat_column} is the sole valid index in the column.")
			print(f"Matched value in the overall flattened and subsetted square is {dummy2dwhere_match}")
		species_config = si.getSpeciesConfig()
		conccount = len(species_config['STATE_VECTOR_CONC'])
		emcount = len(species_config['CONTROL_VECTOR_EMIS'])
		ind_collector = []
		cur_offset = 0
		for i in range(conccount):
			ind_collector.append((dummywhere_match+cur_offset))
			cur_offset+=len(dummywhere_flat)
		for i in range(emcount):
			ind_collector.append((dummy2dwhere_match+cur_offset))
			cur_offset+=len(dummy2dwhere_flat) #Only one value here.
		localizedstatevecinds = np.concatenate(ind_collector)
		if self.verbose>=3:
			print(f"There are a total of {len(localizedstatevecinds)}/{len(self.statevec)} selected from total statevec.")
		return localizedstatevecinds
	def getStateVector(self,latind=None,lonind=None):
		if self.statevec is None:
			self.buildStateVector()
		if not (latind is None): #User supplied ind
			statevecinds = self.getLocalizedStateVectorIndices(latind,lonind)
			statevec_toreturn = self.statevec[statevecinds]
		else: #Return the whole vector
			statevec_toreturn = self.statevec
		if self.verbose>=3:
			print(f"GC Translator number {self.num} got statevector for inds {(latind,lonind)}; this vec has length {len(statevec_toreturn)} of total statevec {len(self.statevec)}.")
		return statevec_toreturn
	#Randomize the restart for purposes of testing. Perturbation is 1/2 of range of percent change selected from a uniform distribution.
	#E.g. 0.1 would range from 90% to 110% of initial values. Bias adds that percent on top of the perturbed fields (0.1 raises everything 10%).
	#Repeats this procedure for every species in the state vector (excluding emissions).
	def randomizeRestart(self,perturbation=0.1,bias=0):
		statevec_species = si.getSpeciesConfig()['STATE_VECTOR_CONC']
		offset = 1-perturbation
		scale = perturbation*2
		for spec in statevec_species:
			conc3d = self.getSpecies3Dconc(spec)
			conc3d *= (scale*np.random.rand(*np.shape(conc3d)))+offset
			conc3d *= 1+bias
			self.setSpecies3Dconc(spec,conc3d)
	#Reconstruct all the 3D concentrations from the analysis vector and overwrite relevant terms in the xr restart dataset.
	#Also construct new scaling factors and add them as a separate array at the new timestep in each of the scaling factor netCDFs.
	#However, only do so for species in the control vectors of emissions and concentrations.
	def reconstructArrays(self,analysis_vector):
		species_config = si.getSpeciesConfig()
		restart_shape = np.shape(self.getSpecies3Dconc(species_config['STATE_VECTOR_CONC'][0]))
		emislist=list(species_config['CONTROL_VECTOR_EMIS'].keys())
		emis_shape = np.shape(self.getEmisSF(emislist[0]))
		counter =  0
		for spec_conc in species_config['STATE_VECTOR_CONC']:
			if spec_conc in species_config['CONTROL_VECTOR_CONC']: #Only overwrite if in the control vector; otherwise just increment.
				index_start = np.sum(self.statevec_lengths[0:counter])
				index_end = np.sum(self.statevec_lengths[0:(counter+1)])
				analysis_subset = analysis_vector[index_start:index_end]
				analysis_3d = np.reshape(analysis_subset,restart_shape) #Unflattens with 'C' order in python
				self.setSpecies3Dconc(spec_conc,analysis_3d) #Overwrite.
			counter+=1
		for spec_emis in species_config['CONTROL_VECTOR_EMIS'].keys(): #Emissions scaling factors are all in the control vector
			index_start = np.sum(self.statevec_lengths[0:counter])
			index_end = np.sum(self.statevec_lengths[0:(counter+1)])
			analysis_subset = analysis_vector[index_start:index_end]
			analysis_emis_2d = np.reshape(analysis_subset,emis_shape) #Unflattens with 'C' order in python
			self.addEmisSF(spec_emis,analysis_emis_2d,species_config['ASSIM_TIME'])
			counter+=1
	def saveRestart(self):
		self.restart_ds["time"] = (["time"], np.array([0]), {"long_name": "Time", "calendar": "gregorian", "axis":"T", "units":self.timestring})
		self.restart_ds.to_netcdf(self.filename)
	def saveEmissions(self):
		for file in self.emis_sf_filenames:
			name = '_'.join(file.split('/')[-1].split('_')[0:-1])
			self.emis_ds_list[name].to_netcdf(file)
