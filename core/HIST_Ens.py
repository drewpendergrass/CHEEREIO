import numpy as np
from glob import glob
import tropomi_tools as tt
import scipy.linalg as la
import toolbox as tx 
import settings_interface as si 
from datetime import date,datetime,timedelta
from HIST_Translator import HIST_Translator

#4D ensemble interface with satellite operators.
class HIST_Ens(object):
	def __init__(self,timestamp,useLevelEdge=False,useStateMet = False,useArea=False,fullperiod=False,interval=None,verbose=1,saveAlbedo=False,useControl=False):
		self.verbose = verbose
		self.saveAlbedo = saveAlbedo
		self.useLevelEdge = useLevelEdge
		self.useStateMet = useStateMet
		self.useArea = useArea
		self.useControl = useControl
		self.spc_config = si.getSpeciesConfig()
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
					self.ht[ens] = HIST_Translator(directory, self.timeperiod,interval,verbose=self.verbose)
				else:
					self.ht[ens] = HIST_Translator(directory, self.timeperiod,verbose=self.verbose)
				ensemble_numbers.append(ens)
		if self.useControl:
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
		self.makeBigY()
	def makeObsTrans(self):
		self.OBS_TRANSLATOR = {}
		self.obsSpecies = []
		for spec in list(self.observed_species.keys()):
			if self.spc_config['OBS_TYPE'][spec]=='TROPOMI':
				self.OBS_TRANSLATOR[spec] = tt.TROPOMI_Translator(self.verbose)
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
	#Gamma^-1 applied in this stage.
	def makeRforSpecies(self,species,latind,lonind):
		errval = float(self.spc_config['OBS_ERROR'][species])
		errtype = self.spc_config['OBS_ERROR_TYPE'][species]
		inds = self.getIndsOfInterest(species,latind,lonind)
		if self.spc_config['AV_TO_GC_GRID']=="True": #If we are averaging to GC grid, the errors will be stored in the ObsData object.
			obsdat = self.bigYDict[species]
			err_av = obsdat.getDataByKey('err_av')
			err_av = err_av[inds]
			to_return = np.diag(err_av**2)
		else:
			if errtype=='absolute':
				to_return = np.diag(np.repeat(errval**2,len(inds))) #we are assuming the user provides the square root of variance
			elif errtype=='relative':
				obsdat = self.bigYDict[species]
				obscol = obsdat.getObsCol()
				obscol = obscol[inds]
				to_return = np.diag((obscol*errval)**2) #multiply measurements by relative error, then square it.
			elif errtype=='product':
				obsdat = self.bigYDict[species]
				err_av = obsdat.getDataByKey('err_av')
				err_av = err_av[inds]
				to_return = np.diag(err_av**2)
		#Apply gamma^-1, so that in the cost function we go from gamma^-1*R to gamma*R^-1
		invgamma = float(self.spc_config['REGULARIZING_FACTOR_GAMMA'][species])**-1
		to_return*=invgamma
		return to_return
	def makeR(self,latind,lonind):
		errmats = []
		for spec in self.obsSpecies:
			errmats.append(self.makeRforSpecies(spec,latind,lonind))
		return la.block_diag(*errmats)
	def getCols(self):
		obsdata_toreturn = {}
		conc2Ds = {}
		errval = float(self.spc_config['OBS_ERROR'][species])
		errcorr = float(self.spc_config['OBS_ERROR_SELF_CORRELATION'][species])
		minerror = float(self.spc_config['MIN_OBS_ERROR'][species])
		errtype = self.spc_config['OBS_ERROR_TYPE'][species]
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
		firstens = self.ensemble_numbers[0]
		hist4D_allspecies = self.ht[firstens].combineHist(self.useLevelEdge,self.useStateMet)
		for species in self.observed_species:
			hist4D = self.ht[firstens].reduceCombinedHistToSpecies(hist4D_allspecies,self.observed_species[species])
			obsdata_toreturn[species] = self.OBS_TRANSLATOR[species].gcCompare(species,self.OBS_DATA[species],hist4D,GC_area=self.AREA,saveAlbedo=self.saveAlbedo,doErrCalc=True,useObserverError=useObserverError,prescribed_error=prescribed_error,prescribed_error_type=prescribed_error_type,transportError = transportError, errorCorr = errcorr,minError=minerror)
			firstcol = obsdata_toreturn[species].getGCCol()
			shape2D = np.zeros(2)
			shape2D[0] = len(firstcol)
			shape2D[1]=len(self.ensemble_numbers)
			shape2D = shape2D.astype(int)
			conc2Ds[species] = np.zeros(shape2D)
			conc2Ds[species][:,firstens-1] = firstcol
		for i in self.ensemble_numbers:
			if i!=firstens:
				hist4D_allspecies = self.ht[i].combineHist(self.useLevelEdge,self.useStateMet)
				for species in self.observed_species:
					hist4D = self.ht[i].reduceCombinedHistToSpecies(hist4D_allspecies,self.observed_species[species])
					col = self.OBS_TRANSLATOR[species].gcCompare(species,self.OBS_DATA[species],hist4D,GC_area=self.AREA,doErrCalc=False).getGCCol()
					conc2Ds[species][:,i-1] = col
		#Save full ensemble data in each of the obsdata objects
		for species in self.observed_species:
			obsdata_toreturn[species].setGCCol(conc2Ds[species])
		if self.useControl:
			hist4D_allspecies = self.control_ht.combineHist(self.useLevelEdge,self.useStateMet)
			for species in self.observed_species:
				hist4D = self.control_ht.reduceCombinedHistToSpecies(hist4D_allspecies,self.observed_species[species])
				col = self.OBS_TRANSLATOR[species].gcCompare(species,self.OBS_DATA[species],hist4D,GC_area=self.AREA,doErrCalc=False).getGCCol()
				obsdata_toreturn[species].addData(control=col)
		return obsdata_toreturn
	def getIndsOfInterest(self,species,latind,lonind):
		loc_rad = float(self.spc_config['LOCALIZATION_RADIUS_km'])
		origlat,origlon = si.getLatLonVals(self.spc_config)
		latval = origlat[latind]
		lonval = origlon[lonind]
		alllat,alllon = self.bigYDict[species].getLatLon()
		distvec = np.array([tx.calcDist_km(latval,lonval,a,b) for a,b in zip(alllat,alllon)])
		inds = np.where(distvec<=loc_rad)[0]
		if len(inds) > self.maxobs:
			inds = np.random.choice(inds, self.maxobs,replace=False) #Randomly subset down to appropriate number of observations
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
