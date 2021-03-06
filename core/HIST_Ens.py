import numpy as np
from glob import glob
import tropomi_tools as tt
import scipy.linalg as la
import toolbox as tx 
from datetime import date,datetime,timedelta
from HIST_Translator import HIST_Translator

#4D ensemble interface with satellite operators.
class HIST_Ens(object):
	def __init__(self,timestamp,useLevelEdge=False,useStateMet = False,useArea=False,fullperiod=False,interval=None,verbose=1,saveAlbedo=False):
		self.verbose = verbose
		self.saveAlbedo = saveAlbedo
		self.useLevelEdge = useLevelEdge
		self.useStateMet = useStateMet
		self.useArea = useArea
		self.spc_config = tx.getSpeciesConfig()
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
				self.SAT_TRANSLATOR[spec] = tt.TROPOMI_Translator(self.verbose)
				self.satSpecies.append(spec)
	def getSatData(self):
		self.SAT_DATA = {}
		for spec in self.satSpecies:
			speciesind = self.satSpecies.index(spec)
			errtype = self.spc_config['OBS_COVARIANCE_TYPE'][speciesind]
			calcError = errtype == 'product'
			self.SAT_DATA[spec] = self.SAT_TRANSLATOR[spec].getObservations(spec,self.timeperiod,self.interval,calcError=calcError)
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
			satdat = self.bigYDict[species]
			satcol = satdat[1]
			satcol = satcol[inds]
			to_return = np.diag((satcol*errval)**2) #multiply measurements by relative error, then square it.
		elif errtype=='product':
			satdat = self.bigYDict[species]
			err_av = satdat[-1] #Will be the last entry
			err_av = err_av[inds]
			to_return = np.diag(err_av**2)
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
		speciesind = self.satSpecies.index(species)
		errval = float(self.spc_config['OBS_COVARIANCE'][speciesind])
		errcorr = float(self.spc_config['OBS_ERROR_SELF_CORRELATION'][speciesind])
		errtype = self.spc_config['OBS_COVARIANCE_TYPE'][speciesind]
		if errtype=='product':
			useError = True
		else:
			useError = False
		firstens = self.ensemble_numbers[0]
		hist4D = self.ht[firstens].combineHist(species,self.useLevelEdge,self.useStateMet)
		if self.spc_config['AV_TO_GC_GRID']=="True":
			if self.saveAlbedo:
				if useError:
					firstcol,satcol,satlat,satlon,sattime,numav,swir_av,nir_av,blended_av,err_av = self.SAT_TRANSLATOR[species].gcCompare(species,self.SAT_DATA[species],hist4D,GC_area=self.AREA,saveAlbedo=True,saveError=True,transportError = errval,errorCorr = errcorr)
				else:
					firstcol,satcol,satlat,satlon,sattime,numav,swir_av,nir_av,blended_av = self.SAT_TRANSLATOR[species].gcCompare(species,self.SAT_DATA[species],hist4D,GC_area=self.AREA,saveAlbedo=True,saveError=False)
			else:
				if useError:
					firstcol,satcol,satlat,satlon,sattime,numav,err_av = self.SAT_TRANSLATOR[species].gcCompare(species,self.SAT_DATA[species],hist4D,GC_area=self.AREA,saveError=True,transportError = errval,errorCorr = errcorr)
				else:
					firstcol,satcol,satlat,satlon,sattime,numav = self.SAT_TRANSLATOR[species].gcCompare(species,self.SAT_DATA[species],hist4D,GC_area=self.AREA,saveError=False)
		else:
			if self.saveAlbedo:
				if useError:
					firstcol,satcol,satlat,satlon,sattime,swir_av,nir_av,blended_av,err_av = self.SAT_TRANSLATOR[species].gcCompare(species,self.SAT_DATA[species],hist4D,GC_area=self.AREA,saveAlbedo=True,saveError=True,transportError = errval,errorCorr = errcorr)
				else:
					firstcol,satcol,satlat,satlon,sattime,swir_av,nir_av,blended_av = self.SAT_TRANSLATOR[species].gcCompare(species,self.SAT_DATA[species],hist4D,GC_area=self.AREA,saveAlbedo=True,saveError=False)
			else:
				if useError:
					firstcol,satcol,satlat,satlon,sattime,err_av = self.SAT_TRANSLATOR[species].gcCompare(species,self.SAT_DATA[species],hist4D,GC_area=self.AREA,saveError=True,transportError = errval,errorCorr = errcorr)
				else:
					firstcol,satcol,satlat,satlon,sattime = self.SAT_TRANSLATOR[species].gcCompare(species,self.SAT_DATA[species],hist4D,GC_area=self.AREA,saveError=False)
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
					col,_,_,_,_,_ = self.SAT_TRANSLATOR[species].gcCompare(species,self.SAT_DATA[species],hist4D,GC_area=self.AREA)
				else:
					col,_,_,_,_ = self.SAT_TRANSLATOR[species].gcCompare(species,self.SAT_DATA[species],hist4D,GC_area=self.AREA)
				conc2D[:,i-1] = col
		if self.spc_config['AV_TO_GC_GRID']=="True":
			if self.saveAlbedo:
				to_return = [conc2D,satcol,satlat,satlon,sattime,numav,swir_av,nir_av,blended_av]
			else:
				to_return = [conc2D,satcol,satlat,satlon,sattime,numav]
		else:
			if self.saveAlbedo:
				to_return = [conc2D,satcol,satlat,satlon,sattime,swir_av,nir_av,blended_av]
			else:
				to_return = [conc2D,satcol,satlat,satlon,sattime]
		if useError:
			to_return.append(err_av)
		return to_return
	def getIndsOfInterest(self,species,latind,lonind):
		loc_rad = float(self.spc_config['LOCALIZATION_RADIUS_km'])
		origlat,origlon = tx.getLatLonVals(self.spc_config)
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
			speciesind = self.satSpecies.index(spec)
			errtype = self.spc_config['OBS_COVARIANCE_TYPE'][speciesind]
			useError = errtype=='product'
			if self.spc_config['AV_TO_GC_GRID']=="True":
				if useError:
					gccol,satcol,_,_,_,_,_ = self.bigYDict[spec]
				else:
					gccol,satcol,_,_,_,_ = self.bigYDict[spec]
			else:
				if useError:
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
