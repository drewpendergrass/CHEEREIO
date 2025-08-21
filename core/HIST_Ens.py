import numpy as np
from glob import glob
import tropomi_tools as tt
import omi_tools as ot
import scipy.linalg as la
from scipy.stats import linregress
import toolbox as tx 
import settings_interface as si
import os.path
from datetime import date,datetime,timedelta
from HIST_Translator import HIST_Translator

#USER: if you have implemented a new observation operator, add it to the operators.json file following the instructions on the Observations page in the documentation
translators = si.importObsTranslators()

#4D ensemble interface with satellite operators.
class HIST_Ens(object):
	def __init__(self,timestamp,useLevelEdge=False,useStateMet = False,useObsPack=False,useArea=False,useSatDiagn=False,fullperiod=False,interval=None,verbose=1,useControl=False):
		self.verbose = verbose
		self.useLevelEdge = useLevelEdge
		self.useStateMet = useStateMet
		self.useObsPack = useObsPack
		self.useArea = useArea
		self.useSatDiagn = useSatDiagn
		self.useControl = useControl
		self.spc_config = si.getSpeciesConfig()
		if self.verbose >=2:
			print(f'HIST_Ens constructor called with the following arguments: HIST_Ens({timestamp},useLevelEdge={useLevelEdge},useStateMet={useStateMet},useObsPack={useObsPack},useArea={useArea},useSatDiagn={useSatDiagn},fullperiod={fullperiod},interval={interval},verbose={verbose},useControl={useControl})')
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
		if self.verbose >=2:
			print(f'HIST_Ens considering observations between {starttime.strftime("%Y/%m/%d %H:%M:%S")} and {endtime.strftime("%Y/%m/%d %H:%M:%S")}')
		self.control_ht = None
		self.ht = {}
		self.observed_species = self.spc_config['OBSERVED_SPECIES']
		self.assimilate_observation = self.spc_config['ASSIMILATE_OBS']
		for ao in self.assimilate_observation:
			self.assimilate_observation[ao] = self.assimilate_observation[ao]=="True" #parse as booleans
		for ens, directory in zip(subdir_numbers,subdirs):
			if (ens==0) and self.useControl:
				if fullperiod:
					self.control_ht = HIST_Translator(directory, self.timeperiod,interval,verbose=self.verbose)
				else:
					self.control_ht = HIST_Translator(directory, self.timeperiod,verbose=self.verbose)
			elif (ens==0) and (not self.useControl): # we don't want to use this directory (because we aren't processing control), but must explicitly omit from processing
				pass #do nothing
			else:
				if fullperiod:
					self.ht[ens] = HIST_Translator(directory, self.timeperiod,interval,verbose=self.verbose)
				else:
					self.ht[ens] = HIST_Translator(directory, self.timeperiod,verbose=self.verbose)
				ensemble_numbers.append(ens)
		if self.useControl and (self.control_ht is None): #Separate control directory
			directory = f"{self.spc_config['MY_PATH']}/{self.spc_config['RUN_NAME']}/control_run/"
			if fullperiod:
				self.control_ht = HIST_Translator(directory, self.timeperiod,interval,verbose=self.verbose)
			else:
				self.control_ht = HIST_Translator(directory, self.timeperiod,verbose=self.verbose)
		self.ensemble_numbers=np.array(ensemble_numbers)
		self.maxobs=int(self.spc_config['MAXNUMOBS'])
		self.interval=interval
		if self.useArea:
			self.AREA = self.ht[1].getArea()
		else:
			self.AREA = None
		if self.verbose >=2:
			print('HIST_Ens construction complete.')
	def makeObsTrans(self):
		self.OBS_TRANSLATOR = {}
		self.obsSpecies = []
		for spec in list(self.observed_species.keys()):
			obstype = self.spc_config['OBS_TYPE'][spec]
			if obstype in translators:
				self.OBS_TRANSLATOR[spec] = translators[obstype](self.verbose)
			else:
				raise ValueError(f'Observer type {obstype} not recongized. Make sure it is listed correctly in operators.json.')
			self.obsSpecies.append(spec)
	def getObsData(self):
		self.OBS_DATA = {}
		for spec in self.obsSpecies:
			errtype = self.spc_config['OBS_ERROR_TYPE'][spec]
			includeObsError = errtype == 'product'
			self.OBS_DATA[spec] = self.OBS_TRANSLATOR[spec].getObservations(spec,self.timeperiod,self.interval,includeObsError=includeObsError)
	def makeBigY(self):
		self.makeObsTrans()
		self.getObsData()
		self.bigYDict = self.getCols()
	#Gamma^-1 applied in this stage. All calculations on a diagonal (vector) until return time
	def makeRforSpecies(self,species,latind,lonind):
		errval = float(self.spc_config['OBS_ERROR'][species])
		errtype = self.spc_config['OBS_ERROR_TYPE'][species]
		inds,distances = self.getIndsOfInterest(species,latind,lonind,return_dist=True)
		if self.spc_config['AV_TO_GC_GRID'][species]=="True": #If we are averaging to GC grid, the errors will be stored in the ObsData object.
			obsdat = self.bigYDict[species]
			err_av = obsdat.getDataByKey('err_av')
			err_av = err_av[inds]
			to_return = err_av**2
		else:
			if errtype=='absolute':
				to_return = np.repeat(errval**2,len(inds)) #we are assuming the user provides the square root of variance
			elif errtype=='relative':
				obsdat = self.bigYDict[species]
				obscol = obsdat.getObsCol()
				obscol = obscol[inds]
				to_return = (obscol*errval)**2#multiply measurements by relative error, then square it.
			elif errtype=='product':
				obsdat = self.bigYDict[species]
				err_av = obsdat.getDataByKey('err_av')
				err_av = err_av[inds]
				to_return = err_av**2
		#Apply gaspari cohn localization.
		if self.spc_config['smooth_localization_with_gaspari_cohn'].lower()=='true':
			loc_rad = float(self.spc_config['LOCALIZATION_RADIUS_km'])
			gaco = tx.make_gaspari_cohn(loc_rad/2)
			weights = gaco(distances) #will be between 0 and 1, shouldn't have anything at zero.
			weights[weights<=0.001] = 0.001 #Set a floor so inverse doesn't explode
			to_return = (weights*to_return**-1)**-1 #Apply gaspari cohn to inverse of R.
		#Apply gamma^-1, so that in the cost function we go from gamma^-1*R to gamma*R^-1
		invgamma = self.getGamma(species)**-1
		to_return*=invgamma
		return np.diag(to_return) #Return as a diagonal matrix.
	def getGamma(self,species):
		diffburnin = self.spc_config['USE_DIFFERENT_GAMMA_FOR_BURN_IN'][species] == "True"
		doburnin = self.spc_config['SIMPLE_SCALE_AT_END_OF_BURN_IN_PERIOD'] == "true"
		if diffburnin and doburnin: #Check that we are (1) doing burnin and (2) adjusting gamma for the given species
			scalingcompleteflag = f"{self.spc_config['MY_PATH']}/{self.spc_config['RUN_NAME']}/scratch/BURN_IN_SCALING_COMPLETE"
			if os.path.isfile(scalingcompleteflag):
				gamma = float(self.spc_config['REGULARIZING_FACTOR_GAMMA'][species]) #Burn in complete, use regular gamma
			else:
				gamma = float(self.spc_config['GAMMA_FOR_BURN_IN'][species]) #In the burn in period, use special gamma.
		else:
			gamma = float(self.spc_config['REGULARIZING_FACTOR_GAMMA'][species])
		return gamma
	def makeR(self,latind,lonind):
		errmats = []
		for spec in self.obsSpecies:
			if self.assimilate_observation[spec]: #If assimilation is turned on, add it to R.
				errmats.append(self.makeRforSpecies(spec,latind,lonind))
		return la.block_diag(*errmats)
	def calcExtrapolationCoefficients(self,species_to_extrapolate):
		gc_version = float(self.spc_config['GC_VERSION'][0:-2])
		if gc_version>=14.1:
			spcconc_name = "SpeciesConcVV"
		else:
			spcconc_name = "SpeciesConc" #Starting in 14.1 we have to specify VV
		extraps = {}
		for i in self.ensemble_numbers: 
			hist4D_allspecies = self.ht[i].combineHist(False,False,False)
			time = hist4D_allspecies.time.values.astype('datetime64[s]').astype('int')
			extraps[i] = {}
			for species in species_to_extrapolate:
				spec4D = hist4D_allspecies[f'{spcconc_name}_{species}'].values
				extraps[i][species] = np.zeros(spec4D.shape[1::])*np.nan
				for j in range(spec4D.shape[1]):
					for k in range(spec4D.shape[2]):
						for l in range(spec4D.shape[3]):
							slope, intercept, r, p, se = linregress(time,spec4D[:,j,k,l])
							extraps[i][species][j,k,l] = slope*60*60 #convert to per hour rather than per second.
		if self.useControl:
			hist4D_allspecies = self.control_ht.combineHist(False,False,False)
			time = hist4D_allspecies.time.values.astype('datetime64[s]').astype('int')
			extraps['control'] = {}
			for species in species_to_extrapolate:
				spec4D = hist4D_allspecies[f'{spcconc_name}_{species}'].values
				extraps['control'][species] = np.zeros(spec4D.shape[1::])*np.nan
				for j in range(spec4D.shape[1]):
					for k in range(spec4D.shape[2]):
						for l in range(spec4D.shape[3]):
							slope, intercept, r, p, se = linregress(time,spec4D[:,j,k,l])
							extraps['control'][species][j,k,l] = slope*60*60 #convert to per hour rather than per second.
		return extraps
	def getCols(self):
		if self.verbose>=2:
			print('HIST_Ens called getCols().')
		to_postprocess = {}
		firstens = self.ensemble_numbers[0]
		hist4D_allspecies = self.ht[firstens].combineHist(self.useLevelEdge,self.useStateMet,self.useObsPack,self.useSatDiagn)
		if self.verbose>=3:
			print('Within getCols(), hist4D produced for the first ensemble member. Details are below:')
			hist4D_allspecies.info()
		for species in self.observed_species:
			to_postprocess[species] = {}
			#Prep all the error information, but just do it once for computational efficiency.
			errval = float(self.spc_config['OBS_ERROR'][species])
			errcorr = float(self.spc_config['OBS_ERROR_SELF_CORRELATION'][species])
			minerror = float(self.spc_config['MIN_OBS_ERROR'][species])
			errtype = self.spc_config['OBS_ERROR_TYPE'][species]
			additional_err_params = self.spc_config['OTHER_OBS_ERROR_PARAMETERS'][species]
			if errtype=='product':
				useObserverError = True
				transportError = errval
				prescribed_error = None
				prescribed_error_type = None
			else:
				useObserverError = False
				transportError = None
				prescribed_error = errval
				prescribed_error_type = errtype
			if "transport_error" in additional_err_params.keys():
				transportError = float(additional_err_params["transport_error"])
			else:
				transportError = None
			gccompare_kwargs = {"GC_area":self.AREA,"doErrCalc":True,"useObserverError":useObserverError,"prescribed_error":prescribed_error,"prescribed_error_type":prescribed_error_type,"transportError":transportError, "errorCorr":errcorr,"minError":minerror}
			#Perform full retrieval with errors, just once
			to_postprocess[species][firstens] = self.OBS_TRANSLATOR[species].gcCompare(species,self.OBS_DATA[species],hist4D_allspecies,**gccompare_kwargs)
			if self.verbose>=3:
				print(f'Within getCols() and for species {species} in the first ensemble member, ObsData generated from Observation Translator gcCompare function with call gcCompare(species={species},OBSDATA=withheld,hist4D_allspecies=withheld,{",".join([f"{key}={gccompare_kwargs[key]}" for key in gccompare_kwargs])})')
		#Perform the rest of the retrievals
		for i in self.ensemble_numbers:
			if i!=firstens:
				hist4D_allspecies = self.ht[i].combineHist(self.useLevelEdge,self.useStateMet,self.useObsPack,self.useSatDiagn)
				if self.verbose>=3:
					print(f'Within getCols(), hist4D produced for ensemble member number {i}. Details are below:')
					hist4D_allspecies.info()
				for species in self.observed_species:
					to_postprocess[species][i] = self.OBS_TRANSLATOR[species].gcCompare(species,self.OBS_DATA[species],hist4D_allspecies,GC_area=self.AREA,doErrCalc=False)
		#Retrieve control, if user asks
		if self.useControl:
			hist4D_allspecies = self.control_ht.combineHist(self.useLevelEdge,self.useStateMet,self.useObsPack,self.useSatDiagn)
			for species in self.observed_species:
				to_postprocess[species]['control'] = self.OBS_TRANSLATOR[species].gcCompare(species,self.OBS_DATA[species],hist4D_allspecies,GC_area=self.AREA,doErrCalc=False)
		#Now apply filters to combine all these obs data objects into single dataset
		obsdata_toreturn = {}
		for species in self.observed_species:
			obsdata_toreturn[species]=applyPostfilter(to_postprocess[species],self.ensemble_numbers)
		return obsdata_toreturn
	def getIndsOfInterest(self,species,latind,lonind,return_dist=False):
		loc_rad = float(self.spc_config['LOCALIZATION_RADIUS_km'])
		origlat,origlon = si.getLatLonVals(self.spc_config)
		latval = origlat[latind]
		lonval = origlon[lonind]
		alllat,alllon = self.bigYDict[species].getLatLon()
		distvec = np.array([tx.calcDist_km(latval,lonval,a,b) for a,b in zip(alllat,alllon)])
		inds = np.where(distvec<=loc_rad)[0]
		if len(inds) > self.maxobs:
			inds = np.random.choice(inds, self.maxobs,replace=False) #Randomly subset down to appropriate number of observations
		if return_dist:
			return [inds,distvec[inds]]
		else:
			return inds
	def getScaling(self,species):
			gccol,obscol = self.bigYDict[species].getCols()
			obsmean = np.mean(gccol,axis=1)
			scaling = np.mean(obscol)/np.mean(obsmean)
			return scaling
	def getLocObsMeanPertDiff(self,latind,lonind):
		obsmeans = []
		obsperts = []
		obsdiffs = []
		for spec in self.obsSpecies:
			if self.assimilate_observation[spec]: #If assimilation is turned on, add it to R.
				ind = self.getIndsOfInterest(spec,latind,lonind)
				errtype = self.spc_config['OBS_ERROR_TYPE'][spec]
				useError = errtype=='product'
				gccol,obscol = self.bigYDict[spec].getCols()
				gccol = gccol[ind,:]
				obscol = obscol[ind]
				obsmean = np.mean(gccol,axis=1)
				obspert = np.zeros(np.shape(gccol))
				for i in range(np.shape(gccol)[1]):
					obspert[:,i]=gccol[:,i]-obsmean
				obsdiff = obscol-obsmean
				obsmeans.append(obsmean)
				obsperts.append(obspert)
				obsdiffs.append(obsdiff)
		full_obsmeans = np.concatenate(obsmeans)
		full_obsperts = np.concatenate(obsperts,axis = 0)
		full_obsdiffs = np.concatenate(obsdiffs)
		return [full_obsmeans,full_obsperts,full_obsdiffs]

#For a list of observation datasets (including control), combine columns into 2D array. 
#If postfilter data is provided, use it to create single consistent array.
#Supply dictionary of form ensnum:obsdata, where control replaces ensnum for control.
def applyPostfilter(dict_of_obsdatas,ensemble_numbers):
	firstens = ensemble_numbers[0]
	toreturn = dict_of_obsdatas[firstens]
	if 'postfilter' in dict_of_obsdatas[firstens].additional_data:
		#Get elements which all obsdatas agree to keep.
		valid_after_postfilter = np.copy(dict_of_obsdatas[firstens].getDataByKey('postfilter'))
		for key in dict_of_obsdatas:
			if key != firstens:
				valid_after_postfilter = np.logical_and(valid_after_postfilter,dict_of_obsdatas[key].getDataByKey('postfilter'))
		#Apply postfilter to everything except gccol, which will come in next step
		toreturn.obscol = toreturn.obscol[valid_after_postfilter]
		toreturn.obslat = toreturn.obslat[valid_after_postfilter]
		toreturn.obslon = toreturn.obslon[valid_after_postfilter]
		toreturn.obstime = toreturn.obstime[valid_after_postfilter]
		for field in toreturn.additional_data: #Apply filter to every other item with right dimensionality
			to_edit = toreturn.getDataByKey(field)
			if len(to_edit)==len(valid_after_postfilter): 
				toreturn.addData(**{field:to_edit[valid_after_postfilter]})
		#Apply postfilter to all the gccols and combine in 2d
		shape2D = np.zeros(2)
		shape2D[0] = len(toreturn.obscol)
		shape2D[1]=len(ensemble_numbers)
		shape2D = shape2D.astype(int)
		conc2D = np.zeros(shape2D)
		for ensnum in ensemble_numbers:
			conc2D[:,ensnum-1] = dict_of_obsdatas[ensnum].getGCCol()[valid_after_postfilter]
		if 'control' in dict_of_obsdatas:
			control = dict_of_obsdatas['control'].getGCCol()[valid_after_postfilter]
		else:
			control = None
	else:
		#Just combine conc2D and get history
		shape2D = np.zeros(2)
		shape2D[0] = len(toreturn.obscol)
		shape2D[1]=len(ensemble_numbers)
		shape2D = shape2D.astype(int)
		conc2D = np.zeros(shape2D)
		for ensnum in ensemble_numbers:
			conc2D[:,ensnum-1] = dict_of_obsdatas[ensnum].getGCCol()
		if 'control' in dict_of_obsdatas:
			control = dict_of_obsdatas['control'].getGCCol()
		else:
			control = None
	#Store data
	toreturn.setGCCol(conc2D)
	if control is not None:
		toreturn.addData(control=control)
	return toreturn
