import numpy as np
import xarray as xr
from glob import glob
import observation_operators as obs
import tropomi_tools as tt
import scipy.linalg as la
import toolbox as tx 
from datetime import date,datetime,timedelta

#This class contains useful methods for getting data from GEOS-Chem restart files and 
#emissions scaling factor netCDFs. After initialization it contains the necessary data
#and can output it in useful ways to other functions in the LETKF procedure.
class GC_Translator(object):
	def __init__(self, path_to_rundir,timestamp,computeStateVec = False,testing=False):
		#self.latinds,self.loninds = tx.getLatLonList(ensnum)
		self.filename = f'{path_to_rundir}GEOSChem.Restart.{timestamp}z.nc4'
		self.timestamp=timestamp
		self.timestring = f'minutes since {timestamp[0:4]}-{timestamp[4:6]}-{timestamp[6:8]} {timestamp[9:11]}:{timestamp[11:13]}:00'
		self.restart_ds = xr.load_dataset(self.filename)
		self.emis_sf_filenames = glob(f'{path_to_rundir}*_SCALEFACTOR.nc')
		self.testing=testing
		if self.testing:
			self.num = path_to_rundir.split('_')[-1][0:4]
			print(f"GC_translator number {self.num} has been called for directory {path_to_rundir} and restart {self.filename}; construction beginning")
		self.emis_ds_list = {}
		for file in self.emis_sf_filenames:
			name = '_'.join(file.split('/')[-1].split('_')[0:-1])
			self.emis_ds_list[name] = xr.load_dataset(file)
			if self.testing:
				print(f"GC_translator number {self.num} has loaded scaling factors for {name}")
		if computeStateVec:
			self.buildStateVector()
		else:
			self.statevec = None
			self.statevec_lengths = None #Until state vector is initialized this variable is None
		if self.testing:
			print(f"GC_Translator number {self.num} construction complete.")
	#Since only one timestamp, returns in format lev,lat,lon
	def getSpecies3Dconc(self, species):
		da = np.array(self.restart_ds[f'SpeciesRst_{species}']).squeeze()
		if self.testing:
			print(f"GC_Translator number {self.num} got 3D conc for species {species} which are of dimension {np.shape(da)}.")
		return da
	def setSpecies3Dconc(self, species, conc3d):
		baseshape = np.shape(conc3d)
		conc4d = conc3d.reshape(np.concatenate([np.array([1]),baseshape]))
		if self.testing:
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
	#We work with the most recent timestamp. Rest are just for archival purposes.
	def getEmisSF(self, species):
		da = self.emis_ds_list[species]['Scalar']
		return np.array(da)[-1,:,:].squeeze()
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
		if tx.getSpeciesConfig(self.testing)['DO_ENS_SPINUP']=='true':
			START_DATE = tx.getSpeciesConfig(self.testing)['ENS_SPINUP_START']
		else:
			START_DATE = tx.getSpeciesConfig(self.testing)['START_DATE']
		orig_timestamp = f'{START_DATE[0:4]}-{START_DATE[4:6]}-{START_DATE[6:8]}' #Start date from  JSON
		END_DATE = tx.getSpeciesConfig(self.testing)['END_DATE']
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
		if self.testing:
			print("*****************************************************************")
			print(f"GC_Translator number {self.num} is starting build of statevector!")
		species_config = tx.getSpeciesConfig(self.testing)
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
		if self.testing:
			print(f"GC_Translator number {self.num} has built statevector; it is of dimension {np.shape(self.statevec)}.")
			print("*****************************************************************")
	def getLocalizedStateVectorIndices(self,latind,lonind):
		surr_latinds, surr_loninds = tx.getIndsOfInterest(latind,lonind,testing=self.testing)
		if self.testing:
			print(f"GC_Translator is getting localized statevec indices surrounding {(latind,lonind)} (lat/lon inds have shapes {np.shape(surr_latinds)}/{np.shape(surr_loninds)}); Lat inds are {surr_latinds} and lon inds are {surr_loninds}.")
		levcount = len(self.getLev())
		latcount = len(self.getLat())
		loncount = len(self.getLon())
		totalcount = levcount*latcount*loncount
		dummy3d = np.arange(0, totalcount).reshape((levcount,latcount,loncount))
		dummywhere_flat = dummy3d[:,surr_latinds,surr_loninds].flatten()
		if self.testing:
			print(f"Within a flattened 3D dummy cube, {len(dummywhere_flat)} entries are valid.")
		dummy2d = np.arange(0, latcount*loncount).reshape((latcount,loncount))
		dummy2dwhere_flat = dummy2d[surr_latinds,surr_loninds].flatten()
		if self.testing:
			print(f"Within a flattened 2D dummy square, {len(dummy2dwhere_flat)} entries are valid.")
		species_config = tx.getSpeciesConfig(self.testing)
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
		if self.testing:
			print(f"There are a total of {len(statevecinds)}/{len(self.statevec)} selected from total statevec.")
		return statevecinds
	def getColumnIndicesFromFullStateVector(self,latind,lonind):
		if self.testing:
			print(f"GC_Translator is getting column statevec indices FOR FULL VECTOR at {(latind,lonind)}.")
		levcount = len(self.getLev())
		latcount = len(self.getLat())
		loncount = len(self.getLon())
		totalcount = levcount*latcount*loncount
		dummy3d = np.arange(0, totalcount).reshape((levcount,latcount,loncount))
		dummywhere_flat = dummy3d[:,latind,lonind].flatten()
		if self.testing:
			print(f"Within a flattened 3D dummy cube, {len(dummywhere_flat)} entries are valid.")
		dummy2d = np.arange(0, latcount*loncount).reshape((latcount,loncount))
		dummy2dwhere_flat = dummy2d[latind,lonind]
		if self.testing:
			print(f"Within a flattened 2D dummy square, {dummy2dwhere_flat} is sole valid entry.")
		species_config = tx.getSpeciesConfig(self.testing)
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
		if self.testing:
			print(f"There are a total of {len(statevecinds)}/{len(self.statevec)} selected from total statevec.")
		return statevecinds
	def getSpeciesConcIndicesInColumn(self,species):
		levcount = len(self.getLev())
		species_config = tx.getSpeciesConfig(self.testing)
		cur_offset = 0
		for ind,spec in enumerate(species_config['STATE_VECTOR_CONC']):
			if species == spec:
				return np.arange(cur_offset,cur_offset+levcount)
			cur_offset+=levcount
		return None #If loop doesn't terminate we did not find the species
	def getSpeciesEmisIndicesInColumn(self,species):
		levcount = len(self.getLev())
		species_config = tx.getSpeciesConfig(self.testing)
		cur_offset = len(species_config['STATE_VECTOR_CONC'])*levcount
		for ind,spec in enumerate(species_config['CONTROL_VECTOR_EMIS']):
			if species == spec:
				return cur_offset
			cur_offset+=1
		return None #If loop doesn't terminate we did not find the species
	def getColumnIndicesFromLocalizedStateVector(self,latind,lonind):
		surr_latinds, surr_loninds = tx.getIndsOfInterest(latind,lonind,testing=self.testing)
		if self.testing:
			print(f"GC_Translator is getting column statevec indices surrounding {(latind,lonind)} (lat/lon inds have shapes {np.shape(surr_latinds)}/{np.shape(surr_loninds)}); Lat inds are {surr_latinds} and lon inds are {surr_loninds}.")
		levcount = len(self.getLev())
		latcount = len(self.getLat())
		loncount = len(self.getLon())
		totalcount = levcount*latcount*loncount
		dummy3d = np.arange(0, totalcount).reshape((levcount,latcount,loncount))
		dummywhere_flat = dummy3d[:,surr_latinds,surr_loninds].flatten()
		dummywhere_flat_column = dummy3d[:,latind,lonind].flatten()
		dummywhere_match = np.where(np.in1d(dummywhere_flat,dummywhere_flat_column))[0]
		if self.testing:
			print(f"Within a flattened 3D dummy cube, {len(dummywhere_flat_column)} entries are valid in the column.")
			print(f"Matched {len(dummywhere_match)} entries in the overall flattened and subsetted column; values are {dummywhere_match}")
		dummy2d = np.arange(0, latcount*loncount).reshape((latcount,loncount))
		dummy2dwhere_flat = dummy2d[surr_latinds,surr_loninds].flatten()
		dummy2dwhere_flat_column = dummy2d[latind,lonind]
		dummy2dwhere_match = np.where(np.in1d(dummy2dwhere_flat,dummy2dwhere_flat_column))[0]
		if self.testing:
			print(f"Within a flattened 2D dummy square, {dummy2dwhere_flat_column} is the sole valid index in the column.")
			print(f"Matched value in the overall flattened and subsetted square is {dummy2dwhere_match}")
		species_config = tx.getSpeciesConfig(self.testing)
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
		if self.testing:
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
		if self.testing:
			print(f"GC Translator number {self.num} got statevector for inds {(latind,lonind)}; this vec has length {len(statevec_toreturn)} of total statevec {len(self.statevec)}.")
		return statevec_toreturn
	#Randomize the restart for purposes of testing. Perturbation is 1/2 of range of percent change selected from a uniform distribution.
	#E.g. 0.1 would range from 90% to 110% of initial values. Bias adds that percent on top of the perturbed fields (0.1 raises everything 10%).
	#Repeats this procedure for every species in the state vector (excluding emissions).
	def randomizeRestart(self,perturbation=0.1,bias=0):
		statevec_species = tx.getSpeciesConfig(self.testing)['STATE_VECTOR_CONC']
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
		species_config = tx.getSpeciesConfig(self.testing)
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

#A class that takes history files and connects them with the main state vector and observation matrices
class HIST_Translator(object):
	def __init__(self, path_to_rundir,timeperiod,interval=None,testing=False):
		self.testing = testing
		self.spc_config = tx.getSpeciesConfig(self.testing)
		self.hist_dir = f'{path_to_rundir}OutputDir'
		self.timeperiod = timeperiod
		self.interval = interval
	def globSubDir(self,timeperiod,useLevelEdge = False, useStateMet = False):
		specconc_list = glob(f'{self.hist_dir}/GEOSChem.SpeciesConc*.nc4')
		specconc_list.sort()
		ts = [datetime.strptime(spc.split('.')[-2][0:13], "%Y%m%d_%H%M") for spc in specconc_list]
		if self.interval:
			specconc_list = [spc for spc,t in zip(specconc_list,ts) if (t>=timeperiod[0]) and (t<timeperiod[1]) and (t.hour % self.interval == 0)]
		else:
			specconc_list = [spc for spc,t in zip(specconc_list,ts) if (t>=timeperiod[0]) and (t<timeperiod[1])]
		if useStateMet:
			met_list = glob(f'{self.hist_dir}/GEOSChem.StateMet*.nc4')
			met_list.sort()
			met_ts = [datetime.strptime(met.split('.')[-2][0:13], "%Y%m%d_%H%M") for met in met_list]
			met_list = [met for met,t in zip(met_list,met_ts) if (t>=timeperiod[0]) and (t<timeperiod[1])]
		if useLevelEdge:
			le_list = glob(f'{self.hist_dir}/GEOSChem.LevelEdgeDiags*.nc4')
			le_list.sort()
			le_ts = [datetime.strptime(le.split('.')[-2][0:13], "%Y%m%d_%H%M") for le in le_list]
			le_list = [le for le,t in zip(le_list,le_ts) if (t>=timeperiod[0]) and (t<timeperiod[1])]
		if useStateMet and useLevelEdge:
			return [specconc_list,le_list,met_list]
		elif useStateMet:
			return [specconc_list,met_list]
		elif useLevelEdge:
			return [specconc_list,le_list]
		else:
			return specconc_list
	def combineHist(self,species,useLevelEdge=False, useStateMet = False):
		dataset=[]
		if useLevelEdge and useStateMet:
			specconc_list,le_list,met_list=self.globSubDir(self.timeperiod,useLevelEdge,useStateMet)
			for specfile,lefile,metfile in zip(specconc_list,le_list,met_list):
				hist_val = xr.load_dataset(specfile)[f'SpeciesConc_{species}']
				lev_val = xr.load_dataset(lefile)[f'Met_PEDGE']
				met_val = xr.load_dataset(metfile)[f'Met_AD']
				data_val = xr.merge([hist_val, lev_val,met_val])
				dataset.append(data_val)
		elif useLevelEdge:
			specconc_list,le_list=self.globSubDir(self.timeperiod,useLevelEdge,useStateMet)
			for specfile,lefile in zip(specconc_list,le_list):
				hist_val = xr.load_dataset(specfile)[f'SpeciesConc_{species}']
				lev_val = xr.load_dataset(lefile)[f'Met_PEDGE']
				data_val = xr.merge([hist_val, lev_val])
				dataset.append(data_val)
		elif useStateMet:
			specconc_list,met_list=self.globSubDir(self.timeperiod,useLevelEdge,useStateMet)
			for specfile,metfile in zip(specconc_list,met_list):
				hist_val = xr.load_dataset(specfile)[f'SpeciesConc_{species}']
				met_val = xr.load_dataset(metfile)[f'Met_AD']
				data_val = xr.merge([hist_val, met_val])
				dataset.append(data_val)
		else:
			specconc_list=self.globSubDir(self.timeperiod,useLevelEdge,useStateMet)
			for specfile in specconc_list:
				hist_val = xr.load_dataset(specfile)[f'SpeciesConc_{species}']
				dataset.append(hist_val)
		dataset = xr.merge(dataset)
		return dataset
	def getArea(self):
		specconc_list=self.globSubDir(self.timeperiod,useLevelEdge=False,useStateMet=False)
		AREA = xr.load_dataset(specconc_list[0])[f'AREA']
		return AREA

#4D ensemble interface with satellite operators.
class HIST_Ens(object):
	def __init__(self,timestamp,useLevelEdge=False,useStateMet = False,useArea=False,fullperiod=False,interval=None,testing=False,saveAlbedo=False):
		self.testing = testing
		self.saveAlbedo = saveAlbedo
		self.useLevelEdge = useLevelEdge
		self.useStateMet = useStateMet
		self.useArea = useArea
		self.spc_config = tx.getSpeciesConfig(self.testing)
		path_to_ensemble = f"{self.spc_config['MY_PATH']}/{self.spc_config['RUN_NAME']}/ensemble_runs"
		subdirs = glob(f"{path_to_ensemble}/*/")
		subdirs.remove(f"{path_to_ensemble}/logs/")
		dirnames = [d.split('/')[-2] for d in subdirs]
		subdir_numbers = [int(n.split('_')[-1]) for n in dirnames]
		ensemble_numbers = []
		endtime = datetime.strptime(timestamp, "%Y%m%d_%H%M")
		if fullperiod:
			START_DATE = self.spc_config['START_DATE']
			starttime = datetime.strptime(f'{START_DATE}_0000', "%Y%m%d_%H%M")
		else:
			ASSIM_TIME = self.spc_config['ASSIM_TIME']
			delta = timedelta(hours=int(ASSIM_TIME))
			starttime = endtime-delta
		self.timeperiod = (starttime,endtime)
		self.ht = {}
		self.observed_species = self.spc_config['OBSERVED_SPECIES']
		for ens, directory in zip(subdir_numbers,subdirs):
			if ens!=0:
				if fullperiod:
					self.ht[ens] = HIST_Translator(directory, self.timeperiod,interval,testing=self.testing)
				else:
					self.ht[ens] = HIST_Translator(directory, self.timeperiod,testing=self.testing)
				ensemble_numbers.append(ens)
		self.ensemble_numbers=np.array(ensemble_numbers)
		self.maxobs=int(self.spc_config['MAXNUMOBS'])
		self.interval=interval
		if self.useArea:
			self.AREA = self.ht[1].getArea()
		else:
			self.AREA = None
		self.makeBigY()
	def makeSatTrans(self):
		self.SAT_TRANSLATOR = {}
		self.satSpecies = []
		for spec,bool4D,boolTROPOMI in zip(list(self.observed_species.values()),self.spc_config['OBS_4D'],self.spc_config['OBS_TYPE_TROPOMI']):
			if (bool4D and boolTROPOMI):
				self.SAT_TRANSLATOR[spec] = tt.TROPOMI_Translator(self.testing)
				self.satSpecies.append(spec)
	def getSatData(self):
		self.SAT_DATA = {}
		for spec in self.satSpecies:
			self.SAT_DATA[spec] = self.SAT_TRANSLATOR[spec].getSatellite(spec,self.timeperiod,self.interval)
	def makeBigY(self):
		self.makeSatTrans()
		self.getSatData()
		self.bigYDict = {}
		for spec in self.satSpecies:
			self.bigYDict[spec] = self.getColsforSpecies(spec)
	#Gamma^-1 applied in this stage.
	def makeRforSpecies(self,species,latind,lonind):
		speciesind = self.satSpecies.index(species)
		errval = float(self.spc_config['OBS_COVARIANCE'][speciesind])
		errtype = self.spc_config['OBS_COVARIANCE_TYPE'][speciesind]
		inds = self.getIndsOfInterest(species,latind,lonind)
		if errtype=='absolute':
			to_return = np.diag(np.repeat(errval,len(inds))) #we are assuming the user already squared these
		elif errtype=='relative':
			if self.spc_config['AV_TO_GC_GRID']=="True":
				_,satcol,_,_,_,_ = self.bigYDict[species]
			else:
				_,satcol,_,_,_ = self.bigYDict[species]
			satcol = satcol[inds]
			to_return = np.diag((satcol*errval)**2) #multiply measurements by relative error, then square it.
		#Apply gamma^-1, so that in the cost function we go from gamma^-1*R to gamma*R^-1
		invgamma = float(self.spc_config['REGULARIZING_FACTOR_GAMMA'][speciesind])**-1
		to_return*=invgamma
		return to_return
	def makeR(self,latind,lonind):
		errmats = []
		for spec in self.satSpecies:
			errmats.append(self.makeRforSpecies(spec,latind,lonind))
		return la.block_diag(*errmats)
	def getColsforSpecies(self,species):
		col3D = []
		firstens = self.ensemble_numbers[0]
		hist4D = self.ht[firstens].combineHist(species,self.useLevelEdge,self.useStateMet)
		if self.spc_config['AV_TO_GC_GRID']=="True":
			if self.saveAlbedo:
				firstcol,satcol,satlat,satlon,sattime,numav,swir_av,nir_av,blended_av = self.SAT_TRANSLATOR[species].gcCompare(species,self.timeperiod,self.SAT_DATA[species],hist4D,GC_area=self.AREA,saveAlbedo=True)
			else:
				firstcol,satcol,satlat,satlon,sattime,numav = self.SAT_TRANSLATOR[species].gcCompare(species,self.timeperiod,self.SAT_DATA[species],hist4D,GC_area=self.AREA)
		else:
			if self.saveAlbedo:
				firstcol,satcol,satlat,satlon,sattime,swir_av,nir_av,blended_av = self.SAT_TRANSLATOR[species].gcCompare(species,self.timeperiod,self.SAT_DATA[species],hist4D,GC_area=self.AREA,saveAlbedo=True)
			else:
				firstcol,satcol,satlat,satlon,sattime = self.SAT_TRANSLATOR[species].gcCompare(species,self.timeperiod,self.SAT_DATA[species],hist4D,GC_area=self.AREA)
		shape2D = np.zeros(2)
		shape2D[0] = len(firstcol)
		shape2D[1]=len(self.ensemble_numbers)
		shape2D = shape2D.astype(int)
		conc2D = np.zeros(shape2D)
		conc2D[:,firstens-1] = firstcol
		for i in self.ensemble_numbers:
			if i!=firstens:
				hist4D = self.ht[i].combineHist(species,self.useLevelEdge,self.useStateMet)
				if self.spc_config['AV_TO_GC_GRID']=="True":
					col,_,_,_,_,_ = self.SAT_TRANSLATOR[species].gcCompare(species,self.timeperiod,self.SAT_DATA[species],hist4D,GC_area=self.AREA)
				else:
					col,_,_,_,_ = self.SAT_TRANSLATOR[species].gcCompare(species,self.timeperiod,self.SAT_DATA[species],hist4D,GC_area=self.AREA)
				conc2D[:,i-1] = col
		if self.spc_config['AV_TO_GC_GRID']=="True":
			if self.saveAlbedo:
				return [conc2D,satcol,satlat,satlon,sattime,numav,swir_av,nir_av,blended_av]
			else:
				return [conc2D,satcol,satlat,satlon,sattime,numav]
		else:
			if self.saveAlbedo:
				return [conc2D,satcol,satlat,satlon,sattime,swir_av,nir_av,blended_av]
			else:
				return [conc2D,satcol,satlat,satlon,sattime]
	def getIndsOfInterest(self,species,latind,lonind):
		loc_rad = float(self.spc_config['LOCALIZATION_RADIUS_km'])
		origlat,origlon = tx.getLatLonVals(self.spc_config,self.testing)
		latval = origlat[latind]
		lonval = origlon[lonind]
		distvec = np.array([tx.calcDist_km(latval,lonval,a,b) for a,b in zip(self.bigYDict[species][2],self.bigYDict[species][3])])
		inds = np.where(distvec<=loc_rad)[0]
		if len(inds) > self.maxobs:
			inds = np.random.choice(inds, self.maxobs,replace=False) #Randomly subset down to appropriate number of observations
		return inds
	def getLocObsMeanPertDiff(self,latind,lonind):
		obsmeans = []
		obsperts = []
		obsdiffs = []
		for spec in self.satSpecies:
			ind = self.getIndsOfInterest(spec,latind,lonind)
			if self.spc_config['AV_TO_GC_GRID']=="True":
				gccol,satcol,_,_,_,_ = self.bigYDict[spec]
			else:
				gccol,satcol,_,_,_ = self.bigYDict[spec]
			gccol = gccol[ind,:]
			satcol = satcol[ind]
			obsmean = np.mean(gccol,axis=1)
			obspert = np.zeros(np.shape(gccol))
			for i in range(np.shape(gccol)[1]):
				obspert[:,i]=gccol[:,i]-obsmean
			obsdiff = satcol-obsmean
			obsmeans.append(obsmean)
			obsperts.append(obspert)
			obsdiffs.append(obsdiff)
		full_obsmeans = np.concatenate(obsmeans)
		full_obsperts = np.concatenate(obsperts,axis = 0)
		full_obsdiffs = np.concatenate(obsdiffs)
		return [full_obsmeans,full_obsperts,full_obsdiffs]

#Lightweight container for GC_Translators; used to combine columns, update restarts, and diff columns.
class GT_Container(object):
	def __init__(self,timestamp,testing=False,constructStateVecs=True):
		self.testing = testing
		spc_config = tx.getSpeciesConfig(self.testing)
		path_to_ensemble = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}/ensemble_runs"
		self.path_to_scratch = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}/scratch"
		npy_column_files = glob(f'{self.path_to_scratch}/**/*.npy',recursive=True)
		npy_col_names = [file.split('/')[-1] for file in npy_column_files]
		npy_columns = [np.load(file) for file in npy_column_files]
		self.columns = dict(zip(npy_col_names,npy_columns))
		subdirs = glob(f"{path_to_ensemble}/*/")
		subdirs.remove(f"{path_to_ensemble}/logs/")
		dirnames = [d.split('/')[-2] for d in subdirs]
		subdir_numbers = [int(n.split('_')[-1]) for n in dirnames]
		ensemble_numbers = []
		self.gt = {}
		self.nature = None
		self.observed_species = spc_config['OBSERVED_SPECIES']
		for ens, directory in zip(subdir_numbers,subdirs):
			if ens==0:
				self.nature = GC_Translator(directory, timestamp, constructStateVecs,self.testing)
			else:
				self.gt[ens] = GC_Translator(directory, timestamp, constructStateVecs,self.testing)
				ensemble_numbers.append(ens)
		self.ensemble_numbers=np.array(ensemble_numbers)
	#Gets saved column and compares to the original files
	def constructColStatevec(self,latind,lonind):
		firstens = self.ensemble_numbers[0]
		col1indvec = self.gt[firstens].getColumnIndicesFromFullStateVector(latind,lonind)
		backgroundEnsemble = np.zeros((len(col1indvec),len(self.ensemble_numbers)))
		backgroundEnsemble[:,firstens-1] = self.gt[firstens].statevec[col1indvec]
		for i in self.ensemble_numbers:
			if i!=firstens:
				colinds = self.gt[i].getColumnIndicesFromFullStateVector(latind,lonind)
				backgroundEnsemble[:,i-1] = self.gt[i].statevec[colinds]
		return backgroundEnsemble
	def diffColumns(self,latind,lonind):
		filenames = list(self.columns.keys())
		substr = f'lat_{latind}_lon_{lonind}.npy'
		search = [i for i in filenames if substr in i]
		saved_col = self.columns[search[0]]
		backgroundEnsemble = self.constructColStatevec(latind,lonind)
		diff = saved_col-backgroundEnsemble
		return [saved_col,backgroundEnsemble,diff]
	def compareSpeciesConc(self,species,latind,lonind):
		firstens = self.ensemble_numbers[0]
		colind = self.gt[firstens].getSpeciesConcIndicesInColumn(species)
		saved_col,backgroundEnsemble,diff = self.diffColumns(latind,lonind)
		saved_col = saved_col[colind,:]
		backgroundEnsemble = backgroundEnsemble[colind,:]
		diff = diff[colind,:]
		col1indvec = self.nature.getColumnIndicesFromFullStateVector(latind,lonind)
		naturecol = self.nature.statevec[col1indvec][colind]
		print(f'*********************************** {species} CONCENTRATION COLUMN AT INDEX {(latind,lonind)} ************************************')
		for i in range(np.shape(saved_col)[1]):
			print(f' ')
			print(f'{species} in ensemble member {i+1} had background concentration of {100*(backgroundEnsemble[:,i]/naturecol)}% nature')
			print(f'{species} in ensemble member {i+1} had analysis concentration of {100*(saved_col[:,i]/naturecol)}% nature')
			print(f'This represents a percent difference of {100*(diff[:,i]/backgroundEnsemble[:,i])}%')
			print(f' ')
	def compareSpeciesEmis(self,species,latind,lonind):
		firstens = self.ensemble_numbers[0]
		colind = self.gt[firstens].getSpeciesEmisIndicesInColumn(species)
		saved_col,backgroundEnsemble,diff = self.diffColumns(latind,lonind)
		saved_col = saved_col[colind,:] #Now will just be a vector of length NumEnsemble
		backgroundEnsemble = backgroundEnsemble[colind,:]
		diff = diff[colind,:]
		col1indvec = self.nature.getColumnIndicesFromFullStateVector(latind,lonind)
		naturecol = self.nature.statevec[col1indvec][colind]
		print(f'*********************************** {species} EMISSIONS SCALING AT INDEX {(latind,lonind)} ************************************')
		for i in range(len(saved_col)):
			print(f' ')
			print(f'{species} in ensemble member {i+1} had background emissions scaling of {100*(backgroundEnsemble[i]/naturecol)}% nature')
			print(f'{species} in ensemble member {i+1} had analysis emissions scaling of {100*(saved_col[i]/naturecol)}% nature')
			print(f'This represents a percent difference of {100*(diff[i]/backgroundEnsemble[i])}%')
			print(f' ')
	def reconstructAnalysisEnsemble(self):
		self.analysisEnsemble = np.zeros((len(self.gt[1].getStateVector()),len(self.ensemble_numbers)))
		for name, cols in zip(self.columns.keys(),self.columns.values()):
			split_name = name.split('_')
			latind = int(split_name[-3])
			lonind = int(split_name[-1].split('.')[0])
			colinds = self.gt[1].getColumnIndicesFromFullStateVector(latind,lonind)
			self.analysisEnsemble[colinds,:] = cols
	def updateRestartsAndScalingFactors(self):
		for i in self.ensemble_numbers:
			self.gt[i].reconstructArrays(self.analysisEnsemble[:,i-1])
	def saveRestartsAndScalingFactors(self):
		for i in self.ensemble_numbers:
			self.gt[i].saveRestart()
			self.gt[i].saveEmissions()

#Contains a dictionary referencing GC_Translators for every run directory.
#In the special case where there is a nature run present (with number 0)
#store the nature run in GC_Translator object nature.
#Also contains an observation operator (pass in the class you would like to use) for each species to assimilate.
#Class contains function to calculate relvant assimilation variables.
#SPECIAL NOTE ON FILES: we will be assuming that geos-chem stopped and left a restart at assimilation time in each run directory.
#That restart will be overwritten in place (name not changed) so next run starts from the assimilation state vector.
#Emissions scaling factors are most recent available (one assimilation timestep ago). New values will be appended to netCDF. 
class Assimilator(object):
	def __init__(self,timestamp,ensnum,corenum,testing=False):
		self.testing = testing
		self.ensnum = ensnum
		self.corenum = corenum
		self.latinds,self.loninds = tx.getLatLonList(ensnum,corenum,self.testing)
		if self.testing:
			print(f"Assimilator has been called for ens {self.ensnum} core {self.corenum}; construction beginning")
			print(f"This core will be handling lat and lon values {[(latval,lonval) for latval,lonval in zip(self.latinds,self.loninds)]}")
		spc_config = tx.getSpeciesConfig(self.testing)
		path_to_ensemble = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}/ensemble_runs"
		self.path_to_scratch = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}/scratch"
		self.path_to_logs = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}/ensemble_runs/logs"
		self.parfilename = f'ens_{ensnum}_core_{corenum}_time_{timestamp}'
		subdirs = glob(f"{path_to_ensemble}/*/")
		subdirs.remove(f"{path_to_ensemble}/logs/")
		dirnames = [d.split('/')[-2] for d in subdirs]
		if self.testing:
			print(f"The following ensemble directories were detected: {dirnames}")
		subdir_numbers = [int(n.split('_')[-1]) for n in dirnames]
		ensemble_numbers = []
		self.nature = None
		self.emcount = len(spc_config['CONTROL_VECTOR_EMIS'])
		self.MINNUMOBS = int(spc_config['MINNUMOBS'])
		self.MinimumScalingFactorAllowed = [float(s) for s in spc_config["MinimumScalingFactorAllowed"]]
		self.MaximumScalingFactorAllowed = [float(s) for s in spc_config["MaximumScalingFactorAllowed"]]
		self.InflateScalingsToXOfPreviousStandardDeviation = [float(s) for s in spc_config["InflateScalingsToXOfPreviousStandardDeviation"]]
		self.MaximumScaleFactorRelativeChangePerAssimilationPeriod=[float(s) for s in spc_config["MaximumScaleFactorRelativeChangePerAssimilationPeriod"]]
		self.AveragePriorAndPosterior = spc_config["AveragePriorAndPosterior"] == "True"
		self.SaveLevelEdgeDiags = spc_config["SaveLevelEdgeDiags"] == "True"
		self.SaveStateMet = spc_config["SaveStateMet"] == "True"
		self.SaveArea = spc_config["SaveArea"] == "True"
		self.SaveDOFS = spc_config["SaveDOFS"] == "True"
		self.DOFS_filter = float(spc_config["DOFS_filter"])
		self.PriorWeightinPriorPosteriorAverage = float(spc_config["PriorWeightinPriorPosteriorAverage"])
		self.gt = {}
		self.observed_species = spc_config['OBSERVED_SPECIES']
		if self.testing:
			print(f"Begin creating GC Translators with state vectors.")
		for ens, directory in zip(subdir_numbers,subdirs):
			if ens==0:
				self.nature = GC_Translator(directory, timestamp, False,self.testing)
			else: 
				self.gt[ens] = GC_Translator(directory, timestamp, True,self.testing)
				ensemble_numbers.append(ens)
		self.ensemble_numbers=np.array(ensemble_numbers)
		if self.testing:
			print(f"GC Translators created. Ensemble number list: {self.ensemble_numbers}")
		self.inflation = float(spc_config['INFLATION_FACTOR'])
		self.histens = HIST_Ens(timestamp,useLevelEdge=self.SaveLevelEdgeDiags,useStateMet = self.SaveStateMet,useArea=self.SaveArea,testing=self.testing)
		if self.testing:
			print(f"Assimilator construction complete")
	def getLat(self):
		return self.gt[1].getLat() #Latitude of first ensemble member, who should always exist
	def getLon(self):
		return self.gt[1].getLon()
	def getLev(self):
		return self.gt[1].getLev()
	def combineEnsemble(self,latind=None,lonind=None):
		if self.testing:
			print(f'combineEnsemble called in Assimilator for lat/lon inds {(latind,lonind)}')
		firstens = self.ensemble_numbers[0]
		firstvec = self.gt[firstens].getStateVector(latind,lonind)
		statevecs = np.zeros((len(firstvec),len(self.ensemble_numbers)))
		statevecs[:,firstens-1] = firstvec
		for i in self.ensemble_numbers:
			if i!=firstens:
				statevecs[:,i-1] = self.gt[i].getStateVector(latind,lonind)
		if self.testing:
			print(f'Ensemble combined in Assimilator for lat/lon inds {(latind,lonind)} and has dimensions {np.shape(statevecs)}.')
		return statevecs
	def ensMeanAndPert(self,latval,lonval):
		if self.testing:
			print(f'ensMeanAndPert called in Assimilator for lat/lon inds {(latval,lonval)}')
		statevecs = self.combineEnsemble(latval,lonval)
		state_mean = np.mean(statevecs,axis = 1) #calculate ensemble mean
		bigX = np.zeros(np.shape(statevecs))
		for i in range(np.shape(bigX)[1]):
			bigX[:,i] = statevecs[:,i]-state_mean
		if self.testing:
			print(f'Ensemble mean at {(latval,lonval)} has dimensions {np.shape(state_mean)} and bigX at at {(latval,lonval)} has dimensions {np.shape(bigX)}.')
		return [state_mean,bigX]
	def combineEnsembleForSpecies(self,species):
		if self.testing:
			print(f'combineEnsembleForSpecies called in Assimilator for species {species}')
		conc3D = []
		firstens = self.ensemble_numbers[0]
		first3D = self.gt[firstens].getSpecies3Dconc(species)
		shape4D = np.zeros(4)
		shape4D[0:3] = np.shape(first3D)
		shape4D[3]=len(self.ensemble_numbers)
		shape4D = shape4D.astype(int)
		conc4D = np.zeros(shape4D)
		conc4D[:,:,:,firstens-1] = first3D
		for i in self.ensemble_numbers:
			if i!=firstens:
				conc4D[:,:,:,i-1] = self.gt[i].getSpecies3Dconc(species)
		return conc4D
	def prepareMeansAndPerts(self,latval,lonval):
		if self.testing:
			print(f'prepareMeansAndPerts called in Assimilator for lat/lon inds {(latval,lonval)}')
		self.ybar_background, self.Ypert_background, self.ydiff = self.histens.getLocObsMeanPertDiff(latval,lonval)
		self.xbar_background, self.Xpert_background = self.ensMeanAndPert(latval,lonval)
		if self.testing:
			print(f'ybar_background for lat/lon inds {(latval,lonval)} has shape {np.shape(self.ybar_background)}.')
			print(f'Ypert_background for lat/lon inds {(latval,lonval)} has shape {np.shape(self.Ypert_background)}.')
			print(f'ydiff for lat/lon inds {(latval,lonval)} has shape {np.shape(self.ydiff)}.')
			print(f'xbar_background for lat/lon inds {(latval,lonval)} has shape {np.shape(self.xbar_background)}.')
			print(f'Xpert_background for lat/lon inds {(latval,lonval)} has shape {np.shape(self.Xpert_background)}.')
	def makeR(self,latind=None,lonind=None):
		if self.testing:
			print(f"Making R for lat/lon inds {(latind,lonind)}.")
		self.R = self.histens.makeR(latind,lonind)
		if self.testing:
			print(f'R for {(latind,lonind)} has dimension {np.shape(self.R)} and value {self.R}')
	def makeC(self):
		self.C = np.transpose(self.Ypert_background) @ la.inv(self.R)
		if self.testing:
			print(f'C made in Assimilator. It has dimension {np.shape(self.C)} and value {self.C}')
	def makePtildeAnalysis(self):
		cyb = self.C @ self.Ypert_background
		k = len(self.ensemble_numbers)
		iden = (k-1)*np.identity(k)/(1+self.inflation)
		self.PtildeAnalysis = la.inv(iden+cyb)
		if self.testing:
			print(f'PtildeAnalysis made in Assimilator. It has dimension {np.shape(self.PtildeAnalysis)} and value {self.PtildeAnalysis}')
	def makeWAnalysis(self):
		k = len(self.ensemble_numbers)
		self.WAnalysis = la.sqrtm((k-1)*self.PtildeAnalysis)
		if self.testing:
			print(f'WAnalysis initialized in Assimilator. It has dimension {np.shape(self.WAnalysis)} and value {self.WAnalysis}')
	def makeWbarAnalysis(self):
		self.WbarAnalysis = self.PtildeAnalysis@self.C@self.ydiff
		if self.testing:
			print(f'WbarAnalysis made in Assimilator. It has dimension {np.shape(self.WbarAnalysis)} and value {self.WbarAnalysis}')
	def adjWAnalysis(self):
		k = len(self.ensemble_numbers)
		for i in range(k):
			self.WAnalysis[:,i]+=self.WbarAnalysis
		if self.testing:
			print(f'WAnalysis adjusted in Assimilator. It has dimension {np.shape(self.WAnalysis)} and value {self.WAnalysis}')
	def makeAnalysisCombinedEnsemble(self):
		self.analysisEnsemble = np.zeros(np.shape(self.Xpert_background))
		k = len(self.ensemble_numbers)
		if np.isnan(self.DOFS_filter):
			for i in range(k):
				self.analysisEnsemble[:,i] = self.Xpert_background.dot(self.WAnalysis[:,i])+self.xbar_background
		else:
			self.analysisPertEnsemble = np.zeros(np.shape(self.Xpert_background))
			for i in range(k):
				self.analysisPertEnsemble[:,i] = self.Xpert_background.dot(self.WAnalysis[:,i])
				self.analysisEnsemble[:,i] = self.analysisPertEnsemble[:,i]+self.xbar_background
		if self.testing:
			print(f'analysisEnsemble made in Assimilator. It has dimension {np.shape(self.analysisEnsemble)} and value {self.analysisEnsemble}')
	def getAnalysisAndBackgroundColumn(self,latval,lonval,doBackground=True, doPerts=False):
		colinds = self.gt[1].getColumnIndicesFromLocalizedStateVector(latval,lonval)
		analysisSubset = self.analysisEnsemble[colinds,:]
		if doPerts:
			analysisPertSubset = self.analysisPertEnsemble[colinds,:]
		if doBackground:
			backgroundSubset = np.zeros(np.shape(self.Xpert_background[colinds,:]))
			k = len(self.ensemble_numbers)
			if doPerts:
				backgroundPertSubset = np.zeros(np.shape(self.Xpert_background[colinds,:]))
				for i in range(k):
					backgroundPertSubset[:,i] = self.Xpert_background[colinds,i]
					backgroundSubset[:,i] = backgroundPertSubset[:,i] +self.xbar_background[colinds]
				return [analysisSubset,backgroundSubset,analysisPertSubset,backgroundPertSubset]
			else:
				for i in range(k):
					backgroundSubset[:,i] = self.Xpert_background[colinds,i]+self.xbar_background[colinds]
				return [analysisSubset,backgroundSubset]
		else:
			return analysisSubset
	def calculateDOFS(self, analysisPert,backgroundPert):
		k = len(self.ensemble_numbers)
		prior_err_cov = (1/(k-1))*backgroundPert@np.transpose(backgroundPert)
		#I checked the math and eqn 13 and eqn 23 in Hunt et al are equivalent bc P tilde is symmetric
		#DCP, 28.03.2022
		posterior_err_cov = (1/(k-1))*analysisPert@np.transpose(analysisPert)
		#DOFS per eqn 11.91 in Brasseur and Jacob is n-tr(Posterior @ Prior -1). Prior isn't full rank so use pseudoinverse.
		traceval = np.trace(posterior_err_cov@np.linalg.pinv(prior_err_cov))
		n_elements = np.shape(posterior_err_cov)[0]
		dofs = n_elements-traceval
		return dofs
	#If for whatever reason we want to skip optimization this can be called.
	def setPosteriorEqualToPrior(self,latval,lonval, returnSubset = True):
		self.analysisEnsemble = np.zeros(np.shape(self.Xpert_background))
		k = len(self.ensemble_numbers)
		for i in range(k):
			self.analysisEnsemble[:,i] = self.Xpert_background[:,i]+self.xbar_background
		if returnSubset:
			analysisSubset = self.getAnalysisAndBackgroundColumn(latval,lonval,doBackground=False)	
			return analysisSubset
	def applyAnalysisCorrections(self,analysisSubset,backgroundSubset):
		#Get scalefactors off the end of statevector
		analysisScalefactor = analysisSubset[(-1*self.emcount)::,:] #This is the column being assimilated, so only one emissions factor per species grouping
		backgroundScalefactor = backgroundSubset[(-1*self.emcount)::,:]
		#Inflate scalings to the X percent of the background standard deviation, per Miyazaki et al 2015
		for i in range(len(self.InflateScalingsToXOfPreviousStandardDeviation)):
			inflator = self.InflateScalingsToXOfPreviousStandardDeviation[i]
			if ~np.isnan(inflator):
				analysis_std = np.std(analysisScalefactor[i,:])
				background_std = np.std(backgroundScalefactor[i,:])
				ratio=analysis_std/background_std
				if ~np.isnan(ratio): #Sometimes background standard deviation is approximately 0.
					if ratio < inflator:
						new_std = inflator*background_std
						meanrebalance = np.mean(analysisScalefactor[i,:])*((new_std/analysis_std)-1)
						analysisScalefactor[i,:] = analysisScalefactor[i,:]*(new_std/analysis_std)-meanrebalance #Scale so sd is new_std and mean is old mean
		#Apply maximum relative change per assimilation period:
		for i in range(len(self.MaximumScaleFactorRelativeChangePerAssimilationPeriod)):
			maxchange=self.MaximumScaleFactorRelativeChangePerAssimilationPeriod[i]
			if ~np.isnan(maxchange):
				relativechanges=(analysisScalefactor[i,:]-backgroundScalefactor[i,:])/backgroundScalefactor[i,:]
				relOverwrite = np.where(np.abs(relativechanges)>maxchange)[0]
				analysisScalefactor[i,relOverwrite] = (1+(np.sign(relativechanges[relOverwrite])*maxchange))*backgroundScalefactor[i,relOverwrite]
		#Set min/max scale factor:
		for i in range(len(self.MinimumScalingFactorAllowed)):
			if ~np.isnan(self.MinimumScalingFactorAllowed[i]):
				minOverwrite = np.where(analysisScalefactor[i,:]<self.MinimumScalingFactorAllowed[i])[0]
				analysisScalefactor[i,minOverwrite] = self.MinimumScalingFactorAllowed[i]
			if ~np.isnan(self.MaximumScalingFactorAllowed[i]):
				maxOverwrite = np.where(analysisScalefactor[i,:]>self.MaximumScalingFactorAllowed[i])[0]
				analysisScalefactor[i,maxOverwrite] = self.MaximumScalingFactorAllowed[i]
		#Done with the scalings
		analysisSubset[(-1*self.emcount)::,:] = analysisScalefactor
		#Now average with prior
		if self.AveragePriorAndPosterior:
			priorweight = self.PriorWeightinPriorPosteriorAverage
			if (priorweight<0) or (priorweight>1):
				raise ValueError('Invalid prior weight; must be between 0 and 1.') 
			posteriorweight = 1-priorweight
			analysisSubset = (backgroundSubset*priorweight)+(analysisSubset*posteriorweight)
		return analysisSubset
	def saveColumn(self,latval,lonval,analysisSubset):
		np.save(f'{self.path_to_scratch}/{str(self.ensnum).zfill(3)}/{str(self.corenum).zfill(3)}/{self.parfilename}_lat_{latval}_lon_{lonval}.npy',analysisSubset)
	def LETKF(self):
		if self.testing:
			print(f"LETKF called! Beginning loop.")
		for latval,lonval in zip(self.latinds,self.loninds):
			if self.testing:
				print(f"Beginning LETKF loop for lat/lon inds {(latval,lonval)}.")
			if self.SaveDOFS:
				latlen = len(self.gt[ens].getLat())
				lonlen = len(self.gt[ens].getLon())
				dofsmat = np.nan*np.zeros((latlen,lonlen))
			self.prepareMeansAndPerts(latval,lonval)
			if len(self.ybar_background)<self.MINNUMOBS:
				#If we don't have enough observations, set posterior equal to prior
				analysisSubset = self.setPosteriorEqualToPrior(latval,lonval,returnSubset=True)
				if self.SaveDOFS:
					dofs = -1
			else:
				self.makeR(latval,lonval)
				self.makeC()
				self.makePtildeAnalysis()
				self.makeWAnalysis()
				self.makeWbarAnalysis()
				self.adjWAnalysis()
				self.makeAnalysisCombinedEnsemble()
				if np.isnan(self.DOFS_filter):
					analysisSubset,backgroundSubset = self.getAnalysisAndBackgroundColumn(latval,lonval,doBackground=True,doPerts=False)
					analysisSubset = self.applyAnalysisCorrections(analysisSubset,backgroundSubset)
				else:
					analysisSubset,backgroundSubset,analysisPertSubset,backgroundPertSubset = self.getAnalysisAndBackgroundColumn(latval,lonval,doBackground=True,doPerts=False)
					dofs = self.calculateDOFS(analysisPertSubset,backgroundPertSubset)
					if dofs >= self.DOFS_filter: #DOFS high enough, proceed with corrections and overwrite
						analysisSubset = self.applyAnalysisCorrections(analysisSubset,backgroundSubset) 
					else: #DOFS too low, not enough information to optimize
						analysisSubset=backgroundSubset #set analysis equal to background
			self.saveColumn(latval,lonval,analysisSubset)
			if self.SaveDOFS:
				dofsmat[latval,lonval] = dofs
		if self.SaveDOFS:
			np.save(f'{self.path_to_logs}/{self.parfilename}_dofsgrid.npy',dofsmat)