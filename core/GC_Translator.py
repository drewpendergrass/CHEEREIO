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
		self.emis_sf_filenames = glob(f'{path_to_rundir}*_SCALEFACTOR.nc')
		self.verbose=verbose
		self.species_config = si.getSpeciesConfig()
		self.useLognormal = self.species_config['lognormalErrors'] == "True"
		self.StateVecType = self.species_config['STATE_VECTOR_CONC_REPRESENTATION'] #How concentrations are stored in state vector: 3D, surface, column_sum, or trop_sum
		if self.verbose>=3:
			self.num = path_to_rundir.split('_')[-1][0:4]
			print(f"GC_translator number {self.num} has been called for directory {path_to_rundir} and restart {self.filename}; construction beginning")
		else:
			self.num=None
		self.data = DataBundle(self.filename,self.emis_sf_filenames,self.species_config,self.timestamp,self.timestamp_as_date,self.useLognormal,self.verbose,self.num)
		if computeStateVec:
			self.statevec = StateVector(StateVecType=self.StateVecType,data=self.data,species_config=self.species_config,emis_sf_filenames=self.emis_sf_filenames,verbose=self.verbose,num=self.num)
		else:
			self.statevec = None
		if self.verbose>=3:
			print(f"GC_Translator number {self.num} construction complete.")
	######    BEGIN FUNCTIONS THAT ALIAS DATA BUNDLE    ########
	def getSpecies3Dconc(self, species): #Since only one timestamp, returns in format lev,lat,lon
		return self.data.getSpecies3Dconc(species)
	def setSpecies3Dconc(self, species, conc3d):
		self.data.setSpecies3Dconc(species, conc3d)
	def getTemp(self):
		return self.data.getTemp()
	def getHeight(self):
		return self.data.getHeight()
	def getTropLev(self):
		return self.data.getTropLev()
	def getPressure(self): #gets pressure at midpoints
		return self.data.getPressure()
	def getLat(self):
		return self.data.getLat()
	def getLon(self):
		return self.data.getLon()
	def getLev(self):
		return self.data.getLev()
	def getRestartTime(self):
		return self.data.getRestartTime()
	def getEmisTime(self):
		return self.data.getEmisTime()
	def getEmisSF(self, species): #Get the emissions from the timestamp nearest to the one supplied by the user.
		return self.data.getEmisSF(species)
	def getEmisLat(self, species):
		return self.data.getEmisLat(species)
	def getEmisLon(self, species):
		return self.data.getEmisLon(species)
	def addEmisSF(self, species, emis2d): #Add 2d emissions scaling factors to the end of the emissions scaling factor
		self.data.addEmisSF(species, emis2d)
	######    END FUNCTIONS THAT ALIAS DATA BUNDLE    ########
	######    BEGIN FUNCTIONS THAT ALIAS STATEVECTOR    ########
	def getLocalizedStateVectorIndices(self,latind,lonind):
		return self.statevec.localizeFromFull(latind,lonind,True)
	def getColumnIndicesFromFullStateVector(self,latind,lonind):
		return self.statevec.localizeFromFull(latind,lonind,False)
	def getColumnIndicesFromLocalizedStateVector(self,latind,lonind):
		return self.statevec.localizeFromFull(latind,lonind,'intersect')
	def getStateVector(self,latind=None,lonind=None):
		if self.statevec is None:
			self.statevec = StateVector(StateVecType=self.StateVecType,data=data,species_config = self.species_config,emis_sf_filenames=self.emis_sf_filenames,verbose=self.verbose,num=self.num)
		return self.statevec.getStateVector(latind,lonind)
	######    END FUNCTIONS THAT ALIAS STATEVECTOR    ########
	def setSpeciesConcByColumn(self,species,column2d,useTrop):
		conc3d = self.data.getSpecies3Dconc(species)
		params = {}
		for param in self.statevec.params_needed:
			params[param] = self.data.getMetByCode(param)
		original_column = self.statevec.StateVecFrom3D(conc3d,**params).reshape(np.shape(column2d))
		ratio = column2d/original_column #values by which we scale
		assim3d = conc3d * ratio #scale by ratio
		if useTrop:
			nlev,nlat,nlon = np.shape(conc3d)
			level = np.arange(0,nlev)
			levcube = np.transpose(np.tile(level,nlat*nlon).reshape((nlat,nlon,nlev)),(2,0,1)) #cube of dim lev,lat,lon with value of level index 
			mask = levcube>=params["trop"] #Make mask where levcube is above tropopause level,
			assim3d[mask] = conc3d[mask] # use mask to select strat/meso values from assimilated 3d matrix, and set them to pre-assimilated value (not updating).		
		self.setSpecies3Dconc(species,assim3d)
	#Randomize the restart for purposes of testing. Perturbation is 1/2 of range of percent change selected from a uniform distribution.
	#E.g. 0.1 would range from 90% to 110% of initial values. Bias adds that percent on top of the perturbed fields (0.1 raises everything 10%).
	#Repeats this procedure for every species in the state vector (excluding emissions).
	def randomizeRestart(self,perturbation=0.1,bias=0):
		statevec_species = self.species_config['STATE_VECTOR_CONC']
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
		counter =  0
		for spec_conc in self.species_config['STATE_VECTOR_CONC']:
			if spec_conc in self.species_config['CONTROL_VECTOR_CONC']: #Only overwrite if in the control vector; otherwise just increment.
				restart_shape = np.shape(self.getSpecies3Dconc(spec_conc))
				index_start = np.sum(self.statevec.statevec_lengths[0:counter])
				index_end = np.sum(self.statevec.statevec_lengths[0:(counter+1)])
				analysis_subset = analysis_vector[index_start:index_end]
				if self.StateVecType == '3D': #case 1: analysis subset is 3D concentrations. Just unfurl and overwrite
					analysis_3d = np.reshape(analysis_subset,restart_shape) #Unflattens with 'C' order in python
					self.setSpecies3Dconc(spec_conc,analysis_3d) #Overwrite.
				else:
					analysis_2d = np.reshape(analysis_subset,restart_shape[1::]) #We are in 2D, so unfurl accordingly
					if self.StateVecType == 'surface':
						self.data.setSpeciesConcByLayer(spec_conc, analysis_2d, layer=0) #overwrite surface layer only.
					elif self.StateVecType == 'trop_sum':
						self.setSpeciesConcByColumn(spec_conc,analysis_2d,useTrop=True) #scale all column values within troposphere evenly using column update.
			counter+=1
		for spec_emis in self.species_config['CONTROL_VECTOR_EMIS'].keys(): #Emissions scaling factors are all in the control vector
			emis_shape = np.shape(self.getEmisSF(spec_emis))
			index_start = np.sum(self.statevec.statevec_lengths[0:counter])
			index_end = np.sum(self.statevec.statevec_lengths[0:(counter+1)])
			analysis_subset = analysis_vector[index_start:index_end]
			analysis_emis_2d = np.reshape(analysis_subset,emis_shape) #Unflattens with 'C' order in python
			self.addEmisSF(spec_emis,analysis_emis_2d)
			counter+=1
	def saveRestart(self):
		self.data.restart_ds["time"] = (["time"], np.array([0]), {"long_name": "Time", "calendar": "gregorian", "axis":"T", "units":self.timestring})
		self.data.restart_ds.to_netcdf(self.filename)
	def saveEmissions(self):
		for file in self.emis_sf_filenames:
			name = '_'.join(file.split('/')[-1].split('_')[0:-1])
			if self.useLognormal:
				self.data.emis_ds_list[name]['Scalar'] = np.exp(self.data.emis_ds_list[name]['Scalar']) #Right before saving, transform back to lognormal space if we are using lognormal errors
			self.data.emis_ds_list[name].to_netcdf(file)

#Handles data getting and setting for emissions and concentrations.
#Class exists to prevent mutual dependencies.
class DataBundle(object):
	def __init__(self,rst_filename,emis_sf_filenames,species_config,timestamp_as_string,timestamp_as_date,useLognormal,verbose,num=None):
		self.restart_ds = xr.open_dataset(rst_filename)
		self.species_config = species_config
		self.verbose = verbose
		self.timestamp = timestamp_as_string
		self.timestamp_as_date = timestamp_as_date
		self.useLognormal = useLognormal
		if verbose >= 3:
			self.num = num
		self.emis_ds_list = {}
		for file in emis_sf_filenames:
			name = '_'.join(file.split('/')[-1].split('_')[0:-1])
			self.emis_ds_list[name] = xr.open_dataset(file)
			if self.useLognormal:
				self.emis_ds_list[name]['Scalar'] = np.log(self.emis_ds_list[name]['Scalar']) #If we are using lognormal errors, convert to gaussian space immediately on import.
				if self.verbose>=3:
					print(f"GC_translator number {self.num} has log transformed scaling factors for {name}")
			if self.verbose>=3:
				print(f"GC_translator number {self.num} has loaded scaling factors for {name}")
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
	def setSpeciesConcByLayer(self, species, conc2d, layer):
		da = self.getSpecies3Dconc(species)
		da[layer,:,:] = conc2d #overwrite layer
		self.setSpecies3Dconc(species,da)
	def getTemp(self):
		return np.array(self.restart_ds['Met_TMPU1']).squeeze() #lev, lat ,lon
	def getHeight(self):
		return np.array(self.restart_ds['Met_BXHEIGHT']).squeeze() #lev, lat, lon
	def getTropLev(self):
		return np.array(self.restart_ds['Met_TropLev']).squeeze() #lat, lon
	def getPressure(self): #gets pressure at midpoints
		psurf = np.array(self.restart_ds['Met_PS1DRY']).squeeze() #lat, lon
		hyam = np.array(self.restart_ds['hyam']) #lev
		hybm = np.array(self.restart_ds['hybm']) #lev
		pmid = hyam +  ( psurf[...,None] * hybm ) #Convert with vector operations to pressure at grid midpoints, lat, lon, lev
		pmid = np.transpose(pmid, (2,0,1)) #reorder to lev, lat, lon as expected
		return pmid
	def getMetByCode(self,code):
		if code == "temp":
			to_return = self.getTemp()
		elif code == "pres":
			to_return = self.getPressure()
		elif code == "height":
			to_return = self.getHeight()
		elif code == "trop":
			to_return = self.getTropLev()
		else:
			raise ValueError(f"Code {code} not recognized.")
		return to_return
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
	def addEmisSF(self, species, emis2d):
		timelist = self.getEmisTime()
		last_time = timelist[-1]
		#new_last_time = last_time+np.timedelta64(assim_time,'h') #Add assim time hours to the last timestamp
		tstr = f'{self.timestamp[0:4]}-{self.timestamp[4:6]}-{self.timestamp[6:8]}T{self.timestamp[9:11]}:{self.timestamp[11:13]}:00.000000000'
		new_last_time = np.datetime64(tstr)
		if self.species_config['DO_ENS_SPINUP']=='true':
			START_DATE = self.species_config['ENS_SPINUP_START']
		else:
			START_DATE = self.species_config['START_DATE']
		orig_timestamp = f'{START_DATE[0:4]}-{START_DATE[4:6]}-{START_DATE[6:8]}' #Start date from  JSON
		END_DATE = self.species_config['END_DATE']
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


class StateVector(object):
	def __init__(self,StateVecType,data,species_config,emis_sf_filenames,verbose,num=None):
		self.StateVecType = StateVecType
		self.StateVecFrom3D = MakeStateVecFrom3D(self.StateVecType)
		if self.StateVecType == "3D":
			self.ConcInterp = "3D" #All 3D concentration values present, interpret state vector appropriately.
			self.params_needed = []
		elif self.StateVecType == "surface":
			self.ConcInterp = "2D" #One layer of concentrations effectively present, interpret state vector appropriately.
			self.params_needed = []
		elif self.StateVecType == "column_sum":
			self.ConcInterp = "2D" #One layer of concentrations effectively present, interpret state vector appropriately.
			self.params_needed = ["temp","pres","height"]
		elif self.StateVecType == "trop_sum":
			self.ConcInterp = "2D" #One layer of concentrations effectively present, interpret state vector appropriately.
			self.params_needed = ["temp","pres","height","trop"]
		else:
			raise ValueError(f"State vector type '{self.StateVecType}' not recognized.")
		self.data = data #this is a databundle
		self.verbose = verbose
		if verbose >= 3:
			self.num = num
		#BUILD STATE VECTOR
		if self.verbose>=3:
			print("*****************************************************************")
			print(f"GC_Translator number {self.num} is starting build of statevector!")
		self.species_config = species_config
		statevec_components = []
		params = {}
		for param in self.params_needed:
			params[param] = self.data.getMetByCode(param)
		for spec_conc in self.species_config['STATE_VECTOR_CONC']:
			statevec_components.append(self.StateVecFrom3D(self.data.getSpecies3Dconc(spec_conc),**params))
		#If no scaling factor files, append 1s because this is a nature directory
		if len(emis_sf_filenames)==0:
			lenones = len(self.data.getLat())*len(self.data.getLon())*len(self.species_config['CONTROL_VECTOR_EMIS'])
			statevec_components.append(np.ones(lenones))
		else:
			for spec_emis in self.species_config['CONTROL_VECTOR_EMIS'].keys():
				statevec_components.append(self.data.getEmisSF(spec_emis).flatten())
		self.statevec_lengths = np.array([len(vec) for vec in statevec_components])
		self.statevec = np.concatenate(statevec_components)
		if self.verbose>=3:
			print(f"GC_Translator number {self.num} has built statevector; it is of dimension {np.shape(self.statevec)}.")
			print("*****************************************************************")
	def makeDummy(self,dimension):
		latcount = len(self.data.getLat())
		loncount = len(self.data.getLon())
		if dimension == "3D":
			levcount = len(self.data.getLev())
			dummy = np.arange(0, levcount*latcount*loncount).reshape((levcount,latcount,loncount))
			incrementor = levcount*latcount*loncount
		if dimension == "2D":
			dummy = np.arange(0, latcount*loncount).reshape((latcount,loncount))
			incrementor = latcount*loncount
		return [dummy,incrementor]
	#getSurroundings: True, False, or intersect (for column within localized statevec)
	def getIndices(self,latind,lonind,getSurroundings):
		dummyConc,conc_incrementor = self.makeDummy(self.ConcInterp) #Concentrations can be 3D or 2D depending on user settings
		dummyEmis,emis_incrementor = self.makeDummy("2D")
		if getSurroundings != False:
			surr_latinds, surr_loninds = tx.getIndsOfInterest(latind,lonind)
			if self.ConcInterp == '3D':
				dummy_conc_flat = dummyConc[:,surr_latinds,surr_loninds].flatten()
			else:
				dummy_conc_flat = dummyConc[surr_latinds,surr_loninds].flatten()
			dummy_emis_flat = dummyEmis[surr_latinds,surr_loninds].flatten()
			if getSurroundings == 'intersect':
				if self.ConcInterp == '3D':
					dummy_conc_flat_column = dummyConc[:,latind,lonind].flatten()
				else:
					dummy_conc_flat_column = dummyConc[latind,lonind]
				conc_index = np.where(np.in1d(dummy_conc_flat,dummy_conc_flat_column))[0]
				conc_incrementor = len(dummy_conc_flat) #for these interesecting cases, we increment by the total number of *localized* entries.
				dummy_emis_flat_column = dummyEmis[latind,lonind]
				emis_index = np.where(np.in1d(dummy_emis_flat,dummy_emis_flat_column))[0]
				emis_incrementor = len(dummy_emis_flat)
			else:
				conc_index = dummy_conc_flat
				emis_index = dummy_emis_flat
		else: #no surroundings
			if self.ConcInterp == '3D':
				conc_index = dummyConc[:,latind,lonind].flatten()
			else:
				conc_index = np.array([dummyConc[latind,lonind]])
			emis_index = np.array([dummyEmis[latind,lonind]])
		return [conc_index,conc_incrementor,emis_index,emis_incrementor]
	#old getLocalizedStateVectorIndices is getSurroundings = True; 
	#old getColumnIndicesFromFullStateVector is getSurroundings = False; 
	#old getColumnIndicesFromLocalizedStateVector is getSurroundings = 'intersect'
	def localizeFromFull(self,latind,lonind,getSurroundings):
		conc_index,conc_incrementor,emis_index,emis_incrementor = self.getIndices(latind,lonind,getSurroundings = getSurroundings)
		conccount = len(self.species_config['STATE_VECTOR_CONC'])
		emcount = len(self.species_config['CONTROL_VECTOR_EMIS'])
		ind_collector = []
		cur_offset = 0
		for i in range(conccount):
			ind_collector.append((conc_index+cur_offset))
			cur_offset+=conc_incrementor
		for i in range(emcount):
			ind_collector.append((emis_index+cur_offset))
			cur_offset+=emis_incrementor
		statevecinds = np.concatenate(ind_collector)
		if self.verbose>=3:
			print(f"There are a total of {len(statevecinds)}/{len(self.statevec)} selected from total statevec.")
		return statevecinds
	def getStateVector(self,latind=None,lonind=None):
		if not (latind is None): #User supplied ind
			statevecinds = self.localizeFromFull(latind,lonind,getSurroundings=True)
			statevec_toreturn = self.statevec[statevecinds]
		else: #Return the whole vector
			statevec_toreturn = self.statevec
		if self.verbose>=3:
			print(f"GC Translator number {self.num} got statevector for inds {(latind,lonind)}; this vec has length {len(statevec_toreturn)} of total statevec {len(self.statevec)}.")
		return statevec_toreturn

def MakeStateVecFrom3D(StateVecType):
	if StateVecType == "3D":
		def StateVecFrom3D(conc3D):
			return conc3D.flatten()
	elif StateVecType == "surface":
		def StateVecFrom3D(conc3D):
			return conc3D[0,:,:].flatten()
	elif StateVecType == "column_sum":
		def StateVecFrom3D(conc3D,temp,pres,height):
			#number density (molecules/cm3), we get there from mol/mol to molec/mol (avogadro) -> molec/J (RT) -> molec/m3 (pressure, in Pa) -> molec/cm3 (10^-6)
			nd3D=(((conc3D*6.0221408e23) / (temp*8.31446261815324))) *pres*1e-6
			# partial coumns (molecules/cm2)
			partial_cols=nd3D*height*1e2 #convert from m to cm
			colsum = np.sum(partial_cols,axis=0) #sum up 
			return colsum.flatten()
	elif StateVecType == "trop_sum":
		def StateVecFrom3D(conc3D,temp,pres,height,trop):
			nlev,nlat,nlon = np.shape(conc3D)
			level = np.arange(0,nlev)
			levcube = np.transpose(np.tile(level,nlat*nlon).reshape((nlat,nlon,nlev)),(2,0,1)) #cube of dim lev,lat,lon with value of level index 
			conc3D[levcube>=trop] = 0 #Make mask where levcube is above tropopause level, use it to select strat/meso values from conc3D, and set them to 0.
			#number density (molecules/cm3), we get there from mol/mol to molec/mol (avogadro) -> molec/J (RT) -> molec/m3 (pressure, in Pa) -> molec/cm3 (10^-6)
			nd3D=(((conc3D*6.0221408e23) / (temp*8.31446261815324))) *pres*1e-6
			# partial coumns (molecules/cm2)
			partial_cols=nd3D*height*1e2 #convert from m to cm
			colsum = np.sum(partial_cols,axis=0) #sum up 
			return colsum.flatten()
	else:
		raise ValueError(f"State vector type '{StateVecType}' not recognized.")
	return StateVecFrom3D



