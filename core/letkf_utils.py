import numpy as np
import xarray as xr
from glob import glob
import observation_operators as obs
import scipy.linalg as la
import toolbox as tx 

def getLETKFConfig(testing=False):
	data = tx.getSpeciesConfig(testing)
	err_config = data['OBS_ERROR_MATRICES']
	if '.npy' in err_config[0]: #Load error matrices from numpy files
		raise NotImplementedError 
	else: #Assume list of strings
		errs = np.array([float(e) for e in err_config])
	#Provide a list of observation operator classes in order of the species to assimilate.
	obs_operator_classes = [getattr(obs, s) for s in data['OBS_OPERATORS']]
	#If you are simulating nature (SIMULATE_NATURE=true in setup_ensemble.sh), provide the nature helper class.
	if data['SIMULATE_NATURE'] == "false":
		raise NotImplementedError #No support for real observations yet!
	else:
		nature_h_functions = [getattr(obs, h) for h in data['NATURE_H_FUNCTIONS']]
	inflation = float(data['INFLATION_FACTOR'])
	return [errs, obs_operator_classes,nature_h_functions,inflation]


#This class contains useful methods for getting data from GEOS-Chem restart files and 
#emissions scaling factor netCDFs. After initialization it contains the necessary data
#and can output it in useful ways to other functions in the LETKF procedure.
class GC_Translator(object):
	def __init__(self, path_to_rundir,timestamp,computeStateVec = False,testing=False):
		#self.latinds,self.loninds = tx.getLatLonList(ensnum)
		self.filename = f'{path_to_rundir}GEOSChem.Restart.{timestamp}z.nc4'
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
		new_last_time = last_time+np.timedelta64(assim_time,'h') #Add assim time hours to the last timestamp
		START_DATE = tx.getSpeciesConfig(self.testing)['START_DATE']
		orig_timestamp = f'{START_DATE[0:4]}-{START_DATE[4:6]}-{START_DATE[6:8]}' #Start date from  JSON
		#Create dataset with this timestep's scaling factors
		ds = xr.Dataset(
			{"Scalar": (("time","lat","lon"), np.expand_dims(emis2d,axis = 0),{"long_name": "Scaling factor", "units":"1"})},
			coords={
				"time": (["time"], np.array([new_last_time]), {"long_name": "time", "calendar": "standard", "units":f"hours since {orig_timestamp} 00:00:00"}),
				"lat": (["lat"], self.getEmisLat(species),{"long_name": "Latitude", "units":"degrees_north"}),
				"lon": (["lon"], self.getEmisLon(species),{"long_name": "Longitude", "units":"degrees_east"})
			},
			attrs={
				"Title":"Auto-generated scaling factors",
				"Conventions":"COARDS",
				"History":"Generated by LETKF"
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
		emis_shape = np.shape(self.getEmisSF(species_config['CONTROL_VECTOR_EMIS'].keys()[0]))
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

#Lightweight container for GC_Translators; used to combine columns, update restarts, and diff columns.
class GT_Container(object):
	def __init__(self,timestamp,testing=False,constructStateVecs=False):
		self.testing = testing
		spc_config = tx.getSpeciesConfig(self.testing)
		path_to_ensemble = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}/ensemble_runs"
		self.path_to_scratch = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}/scratch"
		npy_column_files = glob(f'{self.path_to_scratch}/*.npy')
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
		self.parfilename = f'ens_{ensnum}_core_{corenum}_time_{timestamp}'
		subdirs = glob(f"{path_to_ensemble}/*/")
		subdirs.remove(f"{path_to_ensemble}/logs/")
		dirnames = [d.split('/')[-2] for d in subdirs]
		if self.testing:
			print(f"The following ensemble directories were detected: {dirnames}")
		subdir_numbers = [int(n.split('_')[-1]) for n in dirnames]
		ensemble_numbers = []
		self.nature = None
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
		error_multipliers_or_matrices, self.ObsOperatorClass_list,nature_h_functions,self.inflation = getLETKFConfig(self.testing)
		if self.nature is None: #For the time being, we must have a nature run.
			raise NotImplementedError
		else:
			self.NatureHelperInstance = obs.NatureHelper(self.nature,self.observed_species,nature_h_functions,error_multipliers_or_matrices,self.testing)
			self.makeObsOps()
		if self.testing:
			print(f"Assimilator construction complete")
	def getLat(self):
		return self.gt[1].getLat() #Latitude of first ensemble member, who should always exist
	def getLon(self):
		return self.gt[1].getLon()
	def getLev(self):
		return self.gt[1].getLev()
	def makeObsOps(self):
		if self.testing:
			print(f'makeObsOps called in Assimilator')
		self.ObsOp = {}
		for i,obs_spec_key in enumerate(self.observed_species.keys()):
			ObsOp_instance = self.NatureHelperInstance.makeObsOp(obs_spec_key,self.ObsOperatorClass_list[i])
			self.ObsOp[obs_spec_key] = ObsOp_instance
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
	def ensObsMeanPertDiff(self,latval,lonval):
		if self.testing:
			print(f'ensObsMeanPertDiff called in Assimilator for lat/lon inds {(latval,lonval)}')
		obsmeans = []
		obsperts = []
		obsdiffs = []
		for obskey,species in zip(list(self.observed_species.keys()),list(self.observed_species.values())):
			obsmean,obspert  = self.ensObsMeanAndPertForSpecies(obskey,species,latval,lonval)
			obsmeans.append(obsmean)
			obsperts.append(obspert)
			obsdiffs.append(self.obsDiffForSpecies(obskey,obsmean,latval,lonval))
		full_obsmeans = np.concatenate(obsmeans)
		full_obsperts = np.concatenate(obsperts,axis = 0)
		full_obsdiffs = np.concatenate(obsdiffs)
		if self.testing:
			print(f'Full ObsMeans at {(latval,lonval)} has dimensions {np.shape(full_obsmeans)}; Full ObsPerts at {(latval,lonval)} has dimensions {np.shape(full_obsperts)}; and Full ObsDiffs at {(latval,lonval)} has dimensions {np.shape(full_obsdiffs)}.')
		return [full_obsmeans,full_obsperts,full_obsdiffs]
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
	def ensObsMeanAndPertForSpecies(self, observation_key,species,latval,lonval):
		if self.testing:
			print(f'ensObsMeanAndPertForSpecies called for keys {observation_key} -> {species} in Assimilator for lat/lon inds {(latval,lonval)}')
		spec_4D = self.combineEnsembleForSpecies(species)
		return self.ObsOp[observation_key].obsMeanAndPert(spec_4D,latval,lonval)
	def obsDiffForSpecies(self,observation_key,ensvec,latval,lonval):
		if self.testing:
			print(f'prepareMeansAndPerts called for {observation_key} in Assimilator for lat/lon inds {(latval,lonval)}')
		return self.ObsOp[observation_key].obsDiff(ensvec,latval,lonval)
	def prepareMeansAndPerts(self,latval,lonval):
		if self.testing:
			print(f'prepareMeansAndPerts called in Assimilator for lat/lon inds {(latval,lonval)}')
		self.ybar_background, self.Ypert_background, self.ydiff = self.ensObsMeanPertDiff(latval,lonval)
		self.xbar_background, self.Xpert_background = self.ensMeanAndPert(latval,lonval)
		if self.testing:
			print(f'ybar_background for lat/lon inds {(latval,lonval)} has shape {np.shape(self.ybar_background)}.')
			print(f'Ypert_background for lat/lon inds {(latval,lonval)} has shape {np.shape(self.Ypert_background)}.')
			print(f'ydiff for lat/lon inds {(latval,lonval)} has shape {np.shape(self.ydiff)}.')
			print(f'xbar_background for lat/lon inds {(latval,lonval)} has shape {np.shape(self.xbar_background)}.')
			print(f'Xpert_background for lat/lon inds {(latval,lonval)} has shape {np.shape(self.Xpert_background)}.')
	def makeR(self,latval=None,lonval=None):
		if self.testing:
			print(f'makeR called in Assimilator for lat/lon inds {(latval,lonval)}')
		self.R = self.NatureHelperInstance.makeR(latval,lonval)
		if self.testing:
			print(f'R for {(latval,lonval)} has dimension {np.shape(self.R)} and value {self.R}')
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
		for i in range(k):
			self.analysisEnsemble[:,i] = self.Xpert_background.dot(self.WAnalysis[:,i])+self.xbar_background
		if self.testing:
			print(f'analysisEnsemble made in Assimilator. It has dimension {np.shape(self.analysisEnsemble)} and value {self.analysisEnsemble}')
	def saveColumn(self,latval,lonval):
		colinds = self.gt[1].getColumnIndicesFromLocalizedStateVector(latval,lonval)
		analysisSubset = self.analysisEnsemble[colinds,:]
		np.save(f'{self.path_to_scratch}/{self.parfilename}_lat_{latval}_lon_{lonval}.npy',analysisSubset)
	def LETKF(self):
		if self.testing:
			print(f"LETKF called! Beginning loop.")
		for latval,lonval in zip(self.latinds,self.loninds):
			if self.testing:
				print(f"Beginning LETKF loop for lat/lon inds {(latval,lonval)}.")
			self.prepareMeansAndPerts(latval,lonval)
			self.makeR(latval,lonval)
			self.makeC()
			self.makePtildeAnalysis()
			self.makeWAnalysis()
			self.makeWbarAnalysis()
			self.adjWAnalysis()
			self.makeAnalysisCombinedEnsemble()
			self.saveColumn(latval,lonval)