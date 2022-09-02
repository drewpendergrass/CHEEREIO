import numpy as np
import xarray as xr
from glob import glob
import tropomi_tools as tt
import scipy.linalg as la
import settings_interface as si 
from datetime import date,datetime,timedelta
from GC_Translator import GC_Translator
from HIST_Ens import HIST_Ens

#Contains a dictionary referencing GC_Translators for every run directory.
#In the special case where there is a nature run present (with number 0)
#store the nature run in GC_Translator object nature.
#Also contains an observation operator (pass in the class you would like to use) for each species to assimilate.
#Class contains function to calculate relvant assimilation variables.
#SPECIAL NOTE ON FILES: we will be assuming that geos-chem stopped and left a restart at assimilation time in each run directory.
#That restart will be overwritten in place (name not changed) so next run starts from the assimilation state vector.
#Emissions scaling factors are most recent available (one assimilation timestep ago). New values will be appended to netCDF. 
class Assimilator(object):
	def __init__(self,timestamp,ensnum,corenum):
		spc_config = si.getSpeciesConfig()
		self.verbose = int(spc_config['verbose'])
		self.ensnum = ensnum
		self.corenum = corenum
		self.latinds,self.loninds = si.getLatLonList(ensnum,corenum)
		if self.verbose>=2:
			print(f"Assimilator has been called for ens {self.ensnum} core {self.corenum}; construction beginning")
			print(f"This core will be handling lat and lon values {[(latval,lonval) for latval,lonval in zip(self.latinds,self.loninds)]}")
		path_to_ensemble = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}/ensemble_runs"
		self.path_to_scratch = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}/scratch"
		self.path_to_logs = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}/ensemble_runs/logs"
		self.parfilename = f'ens_{ensnum}_core_{corenum}_time_{timestamp}'
		subdirs = glob(f"{path_to_ensemble}/*/")
		subdirs.remove(f"{path_to_ensemble}/logs/")
		dirnames = [d.split('/')[-2] for d in subdirs]
		if self.verbose>=2:
			print(f"The following ensemble directories were detected: {dirnames}")
		subdir_numbers = [int(n.split('_')[-1]) for n in dirnames]
		ensemble_numbers = []
		self.nature = None
		self.emcount = len(spc_config['CONTROL_VECTOR_EMIS'])
		self.MINNUMOBS = int(spc_config['MINNUMOBS'])
		self.MinimumScalingFactorAllowed = spc_config["MinimumScalingFactorAllowed"]
		self.MaximumScalingFactorAllowed = spc_config["MaximumScalingFactorAllowed"]
		self.InflateScalingsToXOfInitialStandardDeviation = spc_config["InflateScalingsToXOfInitialStandardDeviation"]
		self.emis_names = list(spc_config['CONTROL_VECTOR_EMIS'].keys())
		self.InitEmisSTD = {}
		for name in self.emis_names:
			stds = np.load(f'{path_to_ensemble}/{name}_SCALEFACTOR_INIT_STD.npy')
			self.InitEmisSTD[name] = stds
		self.MaximumScaleFactorRelativeChangePerAssimilationPeriod=spc_config["MaximumScaleFactorRelativeChangePerAssimilationPeriod"]
		self.AveragePriorAndPosterior = spc_config["AveragePriorAndPosterior"] == "True"
		self.SaveLevelEdgeDiags = spc_config["SaveLevelEdgeDiags"] == "True"
		self.SaveStateMet = spc_config["SaveStateMet"] == "True"
		self.SaveArea = spc_config["SaveArea"] == "True"
		self.SaveDOFS = spc_config["SaveDOFS"] == "True"
		self.DOFS_filter = float(spc_config["DOFS_filter"])
		self.PriorWeightinPriorPosteriorAverage = float(spc_config["PriorWeightinPriorPosteriorAverage"])
		self.gt = {}
		self.observed_species = spc_config['OBSERVED_SPECIES']
		if self.verbose>=2:
			print(f"Begin creating GC Translators with state vectors.")
		for ens, directory in zip(subdir_numbers,subdirs):
			if ens==0:
				self.nature = GC_Translator(directory, timestamp, False,self.verbose)
			else: 
				self.gt[ens] = GC_Translator(directory, timestamp, True,self.verbose)
				ensemble_numbers.append(ens)
		self.ensemble_numbers=np.array(ensemble_numbers)
		if self.verbose>=2:
			print(f"GC Translators created. Ensemble number list: {self.ensemble_numbers}")
		self.inflation = float(spc_config['INFLATION_FACTOR'])
		self.histens = HIST_Ens(timestamp,useLevelEdge=self.SaveLevelEdgeDiags,useStateMet = self.SaveStateMet,useArea=self.SaveArea,verbose=self.verbose)
		if self.verbose>=2:
			print(f"Assimilator construction complete")
	def getLat(self):
		return self.gt[1].getLat() #Latitude of first ensemble member, who should always exist
	def getLon(self):
		return self.gt[1].getLon()
	def getLev(self):
		return self.gt[1].getLev()
	def combineEnsemble(self,latind=None,lonind=None):
		if self.verbose>=2:
			print(f'combineEnsemble called in Assimilator for lat/lon inds {(latind,lonind)}')
		firstens = self.ensemble_numbers[0]
		firstvec = self.gt[firstens].getStateVector(latind,lonind)
		statevecs = np.zeros((len(firstvec),len(self.ensemble_numbers)))
		statevecs[:,firstens-1] = firstvec
		for i in self.ensemble_numbers:
			if i!=firstens:
				statevecs[:,i-1] = self.gt[i].getStateVector(latind,lonind)
		if self.verbose>=2:
			print(f'Ensemble combined in Assimilator for lat/lon inds {(latind,lonind)} and has dimensions {np.shape(statevecs)}.')
		return statevecs
	def ensMeanAndPert(self,latval,lonval):
		if self.verbose>=2:
			print(f'ensMeanAndPert called in Assimilator for lat/lon inds {(latval,lonval)}')
		statevecs = self.combineEnsemble(latval,lonval)
		state_mean = np.mean(statevecs,axis = 1) #calculate ensemble mean
		bigX = np.zeros(np.shape(statevecs))
		for i in range(np.shape(bigX)[1]):
			bigX[:,i] = statevecs[:,i]-state_mean
		if self.verbose>=2:
			print(f'Ensemble mean at {(latval,lonval)} has dimensions {np.shape(state_mean)} and bigX at at {(latval,lonval)} has dimensions {np.shape(bigX)}.')
		return [state_mean,bigX]
	def combineEnsembleForSpecies(self,species):
		if self.verbose>=2:
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
		if self.verbose>=2:
			print(f'prepareMeansAndPerts called in Assimilator for lat/lon inds {(latval,lonval)}')
		self.ybar_background, self.Ypert_background, self.ydiff = self.histens.getLocObsMeanPertDiff(latval,lonval)
		self.xbar_background, self.Xpert_background = self.ensMeanAndPert(latval,lonval)
		if self.verbose>=2:
			print(f'ybar_background for lat/lon inds {(latval,lonval)} has shape {np.shape(self.ybar_background)}.')
			print(f'Ypert_background for lat/lon inds {(latval,lonval)} has shape {np.shape(self.Ypert_background)}.')
			print(f'ydiff for lat/lon inds {(latval,lonval)} has shape {np.shape(self.ydiff)}.')
			print(f'xbar_background for lat/lon inds {(latval,lonval)} has shape {np.shape(self.xbar_background)}.')
			print(f'Xpert_background for lat/lon inds {(latval,lonval)} has shape {np.shape(self.Xpert_background)}.')
	def makeR(self,latind=None,lonind=None):
		if self.verbose>=2:
			print(f"Making R for lat/lon inds {(latind,lonind)}.")
		self.R = self.histens.makeR(latind,lonind)
		if self.verbose>=2:
			print(f'R for {(latind,lonind)} has dimension {np.shape(self.R)} and value {self.R}')
	def makeC(self):
		self.C = np.transpose(self.Ypert_background) @ la.inv(self.R)
		if self.verbose>=3:
			print(f'C made in Assimilator. It has dimension {np.shape(self.C)} and value {self.C}')
	def makePtildeAnalysis(self):
		cyb = self.C @ self.Ypert_background
		k = len(self.ensemble_numbers)
		iden = (k-1)*np.identity(k)/(1+self.inflation)
		self.PtildeAnalysis = la.inv(iden+cyb)
		if self.verbose>=3:
			print(f'PtildeAnalysis made in Assimilator. It has dimension {np.shape(self.PtildeAnalysis)} and value {self.PtildeAnalysis}')
	def makeWAnalysis(self):
		k = len(self.ensemble_numbers)
		self.WAnalysis = la.sqrtm((k-1)*self.PtildeAnalysis)
		if self.verbose>=3:
			print(f'WAnalysis initialized in Assimilator. It has dimension {np.shape(self.WAnalysis)} and value {self.WAnalysis}')
	def makeWbarAnalysis(self):
		self.WbarAnalysis = self.PtildeAnalysis@self.C@self.ydiff
		if self.verbose>=3:
			print(f'WbarAnalysis made in Assimilator. It has dimension {np.shape(self.WbarAnalysis)} and value {self.WbarAnalysis}')
	def adjWAnalysis(self):
		k = len(self.ensemble_numbers)
		for i in range(k):
			self.WAnalysis[:,i]+=self.WbarAnalysis
		if self.verbose>=3:
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
		if self.verbose>=2:
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
	def applyAnalysisCorrections(self,analysisSubset,backgroundSubset,latind,lonind):
		if self.verbose>=2:
			print(f"applyAnalysisCorrections called for index {(latind,lonind)}.")
		#Get scalefactors off the end of statevector
		analysisScalefactor = analysisSubset[(-1*self.emcount)::,:] #This is the column being assimilated, so only one emissions factor per species grouping
		backgroundScalefactor = backgroundSubset[(-1*self.emcount)::,:]
		if self.verbose>=2:
			print(f"Analysis scale factor has dimension {np.shape(analysisScalefactor)} and value {analysisScalefactor}.")
			print(f"The analysis ensemble mean is {np.mean(analysisScalefactor,axis=1)}.")
			print(f"Background scale factor has dimension {np.shape(backgroundScalefactor)} and value {backgroundScalefactor}.")
			print(f"The background ensemble mean is {np.mean(backgroundScalefactor,axis=1)}.")
		#Inflate scalings to the X percent of the background standard deviation, per Miyazaki et al 2015
		if self.verbose>=2:
			print('BEGIN section InflateScalingsToXOfInitialStandardDeviation')
		for i,emis in enumerate(self.emis_names):
			inflator = float(self.InflateScalingsToXOfInitialStandardDeviation[emis])
			if ~np.isnan(inflator):
				analysis_std = np.std(analysisScalefactor[i,:])
				background_std = self.InitEmisSTD[emis][latind,lonind]
				if self.verbose>=2:
					print(f'Inflating {emis} scalings to {inflator} of initial scaling factor standard deviation.')
					print(f'The {emis} analysis scaling standard deviation is {analysis_std}, against an initial standard deviation of {background_std}.')
				if background_std == 0:
					ratio = np.nan 
				else:
					ratio=analysis_std/background_std
				if ~np.isnan(ratio):
					if self.verbose>=2:
						print(f'This yields a ratio of {ratio}, compared against inflator standard of {inflator}')
					if ratio < inflator:
						new_std = inflator*background_std
						if self.verbose>=2:
							print(f'Ratio {ratio} is less than inflator standard {inflator}. Adjusting data so new standard deviation will be {new_std}.')
						meanrebalance = np.mean(analysisScalefactor[i,:])*((new_std/analysis_std)-1)
						if self.verbose>=2:
							print(f'Old analysis scale factors were {analysisScalefactor[i,:]}, which have mean {np.mean(analysisScalefactor[i,:])} and standard deviation {np.std(analysisScalefactor[i,:])}.')
						analysisScalefactor[i,:] = analysisScalefactor[i,:]*(new_std/analysis_std)-meanrebalance #Scale so sd is new_std and mean is old mean
						if self.verbose>=2:
							print(f'New {emis} analysis scale factors are {analysisScalefactor[i,:]}, which have mean {np.mean(analysisScalefactor[i,:])} and standard deviation {np.std(analysisScalefactor[i,:])}.')
					else:
						if self.verbose>=2:
							print(f'Ratio {ratio} is greater than inflator standard {inflator}, so doing nothing.')
		if self.verbose>=2:
			print('END section InflateScalingsToXOfInitialStandardDeviation')
		if self.verbose>=2:
			print('BEGIN section MaximumScaleFactorRelativeChangePerAssimilationPeriod')
		#Apply maximum relative change per assimilation period:
		for i,emis in enumerate(self.emis_names):
			maxchange = float(self.MaximumScaleFactorRelativeChangePerAssimilationPeriod[emis])
			if ~np.isnan(maxchange):
				if self.verbose>=2:
					print(f'If changes in scaling factor for {emis} is greater than {maxchange}, saturate.')
				relativechanges=(analysisScalefactor[i,:]-backgroundScalefactor[i,:])/backgroundScalefactor[i,:]
				if self.verbose>=2:
					print(f'Relative changes for {emis} are {relativechanges}.')
				relOverwrite = np.where(np.abs(relativechanges)>maxchange)[0]
				if self.verbose>=2:
					print(f'The following locations exceed the saturation cap: {relOverwrite}')
				analysisScalefactor[i,relOverwrite] = (1+(np.sign(relativechanges[relOverwrite])*maxchange))*backgroundScalefactor[i,relOverwrite]
				if self.verbose>=2:
					print(f'New {emis} analysis scale factors are {analysisScalefactor[i,:]}.')
		if self.verbose>=2:
			print('END section MaximumScaleFactorRelativeChangePerAssimilationPeriod')
		if self.verbose>=2:
			print('BEGIN section Minimum/MaximumScalingFactorAllowed')
		#Set min/max scale factor:
		for i,emis in enumerate(self.emis_names):
			minsf = float(self.MinimumScalingFactorAllowed[emis])
			maxsf = float(self.MaximumScalingFactorAllowed[emis])
			if ~np.isnan(minsf):
				minOverwrite = np.where(analysisScalefactor[i,:]<minsf)[0]
				if self.verbose>=2:
					print(f'The following {emis} scaling factors fall below minimum of {minsf}: {analysisScalefactor[i,minOverwrite]}')
				analysisScalefactor[i,minOverwrite] = minsf
			if ~np.isnan(maxsf):
				maxOverwrite = np.where(analysisScalefactor[i,:]>maxsf)[0]
				if self.verbose>=2:
					print(f'The following {emis} scaling factors rise above maximum of {maxsf}: {analysisScalefactor[i,maxOverwrite]}')
				analysisScalefactor[i,maxOverwrite] = maxsf
			if self.verbose>=2:
				print(f'New analysis scale factors are {analysisScalefactor[i,:]}.')
		if self.verbose>=2:
			print('END section Minimum/MaximumScalingFactorAllowed')
		#Done with the scalings
		if self.verbose>=2:
			print(f'Old scaling factors at end of analysis subset: {analysisSubset[(-1*self.emcount)::,:]}')
		analysisSubset[(-1*self.emcount)::,:] = analysisScalefactor
		if self.verbose>=2:
			print(f'New scaling factors at end of analysis subset: {analysisSubset[(-1*self.emcount)::,:]}')
		#Now average with prior
		if self.AveragePriorAndPosterior:
			priorweight = self.PriorWeightinPriorPosteriorAverage
			if (priorweight<0) or (priorweight>1):
				raise ValueError('Invalid prior weight; must be between 0 and 1.') 
			posteriorweight = 1-priorweight
			if self.verbose>=2:
				print(f'Averaging prior and posterior, with prior weight of {priorweight} and posterior weight of {posteriorweight}.')
			analysisSubset = (backgroundSubset*priorweight)+(analysisSubset*posteriorweight)
		return analysisSubset
	def saveColumn(self,latval,lonval,analysisSubset):
		np.save(f'{self.path_to_scratch}/{str(self.ensnum).zfill(3)}/{str(self.corenum).zfill(3)}/{self.parfilename}_lat_{latval}_lon_{lonval}.npy',analysisSubset)
	def scaleRestarts(self):
		if self.verbose>=1:
			print(f"Scaling all restarts to match observations.")
		scale_factors_by_species_key = {}
		for species_key in self.observed_species:
			scale_factors_by_species_key[species_key] = self.histens.getScaling(species_key) #Get the scaling factor to make GC ens mean match obs mean.
		scale_factors_by_species = {} #If we have multiple scale factors for one species (e.g. surface and satellite observations, average the scalings)
		for species_key in self.observed_species:
			species_value = self.observed_species[species_key]
			if species_value in list(scale_factors_by_species.keys())
				scale_factors_by_species[species_value] = np.mean([scale_factors_by_species_key[species_key],scale_factors_by_species[species_value]])
			else:
				scale_factors_by_species[species_value] = scale_factors_by_species_key[species_key]
		for species in scale_factors_by_species:
			scaling_factor = scale_factors_by_species[species]
			for i in self.ensemble_numbers:
				scaled_species = self.gt[i].getSpecies3Dconc(species)*scaling_factor
				self.gt[i].setSpecies3Dconc(species, scaled_species)
		#Save out the resarts
		for i in self.ensemble_numbers:
			self.gt[i].saveRestart()
	def LETKF(self):
		if self.verbose>=2:
			print(f"LETKF called! Beginning loop.")
		if self.SaveDOFS:
			latlen = len(self.gt[1].getLat())
			lonlen = len(self.gt[1].getLon())
			dofsmat = np.nan*np.zeros((latlen,lonlen))
		for latval,lonval in zip(self.latinds,self.loninds):
			if self.verbose>=1:
				print(f"Beginning LETKF loop for lat/lon inds {(latval,lonval)}.")
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
					analysisSubset = self.applyAnalysisCorrections(analysisSubset,backgroundSubset,latval,lonval)
				else:
					analysisSubset,backgroundSubset,analysisPertSubset,backgroundPertSubset = self.getAnalysisAndBackgroundColumn(latval,lonval,doBackground=True,doPerts=True)
					dofs = self.calculateDOFS(analysisPertSubset,backgroundPertSubset)
					if dofs >= self.DOFS_filter: #DOFS high enough, proceed with corrections and overwrite
						analysisSubset = self.applyAnalysisCorrections(analysisSubset,backgroundSubset,latval,lonval) 
					else: #DOFS too low, not enough information to optimize
						analysisSubset=backgroundSubset #set analysis equal to background
			self.saveColumn(latval,lonval,analysisSubset)
			if self.SaveDOFS:
				dofsmat[latval,lonval] = dofs
		if self.SaveDOFS:
			np.save(f'{self.path_to_logs}/dofs_scratch/{self.parfilename}_dofsgrid.npy',dofsmat)