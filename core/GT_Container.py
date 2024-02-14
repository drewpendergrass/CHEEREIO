import numpy as np
from glob import glob
import settings_interface as si 
from datetime import date,datetime,timedelta
from GC_Translator import GC_Translator

#Lightweight container for GC_Translators; used to combine columns, update restarts, and diff columns.
class GT_Container(object):
	def __init__(self,timestamp,getAssimColumns=True, constructStateVecs=True):
		spc_config = si.getSpeciesConfig()
		path_to_ensemble = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}/ensemble_runs"
		self.path_to_scratch = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}/scratch"
		if getAssimColumns:
			npy_column_files = glob(f'{self.path_to_scratch}/**/*.npy',recursive=True)
			npy_col_names = [file.split('/')[-1] for file in npy_column_files]
			npy_columns = [np.load(file) for file in npy_column_files]
			self.columns = dict(zip(npy_col_names,npy_columns))
		else:
			self.columns = None
		subdirs = glob(f"{path_to_ensemble}/*/")
		subdirs.remove(f"{path_to_ensemble}/logs/")
		dirnames = [d.split('/')[-2] for d in subdirs]
		subdir_numbers = [int(n.split('_')[-1]) for n in dirnames]
		ensemble_numbers = []
		self.gt = {}
		self.observed_species = spc_config['OBSERVED_SPECIES']
		for ens, directory in zip(subdir_numbers,subdirs):
			if ens!=0:
				self.gt[ens] = GC_Translator(directory, timestamp, constructStateVecs)
				ensemble_numbers.append(ens)
		self.ensemble_numbers=np.array(ensemble_numbers)
		self.verbose = int(spc_config['verbose'])
		#Check if we need to do bonus inflation, not handled by Assimilator
		self.species_not_in_statevec_to_RTPS = []
		for species in spc_config['species_not_in_statevec_to_RTPS']:
			#Don't inflate species in state_vector_conc.
			if species not in spc_config['STATE_VECTOR_CONC']:
				self.species_not_in_statevec_to_RTPS.append(species)
		if len(self.species_not_in_statevec_to_RTPS)>0:
			self.RTPS_parameter = float(spc_config["RTPS_parameter"])
	#Gets saved column and compares to the original files
	def constructColStatevec(self,latind,lonind):
		firstens = self.ensemble_numbers[0]
		col1indvec = self.gt[firstens].getColumnIndicesFromFullStateVector(latind,lonind)
		backgroundEnsemble = np.zeros((len(col1indvec),len(self.ensemble_numbers)))
		backgroundEnsemble[:,firstens-1] = self.gt[firstens].getStateVector()[col1indvec]
		for i in self.ensemble_numbers:
			if i!=firstens:
				colinds = self.gt[i].getColumnIndicesFromFullStateVector(latind,lonind)
				backgroundEnsemble[:,i-1] = self.gt[i].getStateVector()[colinds]
		return backgroundEnsemble
	def diffColumns(self,latind,lonind):
		filenames = list(self.columns.keys())
		substr = f'lat_{latind}_lon_{lonind}.npy'
		search = [i for i in filenames if substr in i]
		saved_col = self.columns[search[0]]
		backgroundEnsemble = self.constructColStatevec(latind,lonind)
		diff = saved_col-backgroundEnsemble
		return [saved_col,backgroundEnsemble,diff]
	#Same func as in Assimilator
	def combineEnsembleForSpecies(self,species):
		if self.verbose>=2:
			print(f'combineEnsembleForSpecies called in Assimilator for species {species}')
		conc3D = []
		firstens = self.ensemble_numbers[0]
		first3D = self.gt[firstens].getSpecies3Dconc(species)
		shape4D = np.zeros(4)
		shape4D[1:4] = np.shape(first3D)
		shape4D[0]=len(self.ensemble_numbers)
		shape4D = shape4D.astype(int)
		conc4D = np.zeros(shape4D)
		conc4D[firstens-1,:,:,:] = first3D
		for i in self.ensemble_numbers:
			if i!=firstens:
				conc4D[i-1,:,:,:] = self.gt[i].getSpecies3Dconc(species)
		return conc4D
	def constructBackgroundEnsemble(self):
		self.backgroundEnsemble = np.zeros((len(self.gt[1].getStateVector()),len(self.ensemble_numbers)))
		for i in self.ensemble_numbers:
			self.backgroundEnsemble[:,i-1] = self.gt[i].getStateVector()
	def reconstructAnalysisEnsemble(self):
		self.analysisEnsemble = np.zeros((len(self.gt[1].getStateVector()),len(self.ensemble_numbers)))
		for name, cols in zip(self.columns.keys(),self.columns.values()):
			split_name = name.split('_')
			latind = int(split_name[-3])
			lonind = int(split_name[-1].split('.')[0])
			colinds = self.gt[1].getColumnIndicesFromFullStateVector(latind,lonind)
			self.analysisEnsemble[colinds,:] = cols
	#For some forms of additional inflation requested by user, we need to store out data from the background state
	def prepForInflation(self):
		#For RTPS, we need to save out background std
		if len(self.species_not_in_statevec_to_RTPS)>0:
			self.sigma_b = {}
			for species in self.species_not_in_statevec_to_RTPS:
				conc4D = self.combineEnsembleForSpecies(species)
				self.sigma_b[species] = np.std(conc4D,axis=0) 
	def updateRestartsAndScalingFactors(self):
		for i in self.ensemble_numbers:
			self.gt[i].reconstructArrays(self.analysisEnsemble[:,i-1])
	#In some cases, user requests inflation performed on species not in state vector. Do that now.
	def performAdditionalInflation(self):
		#Do RTPS for select species not in statevector
		if len(self.species_not_in_statevec_to_RTPS)>0:
			for species in self.species_not_in_statevec_to_RTPS:
				conc4D = self.combineEnsembleForSpecies(species)
				sigma_a = np.std(conc4D,axis=0)
				sigma_RTPS = (self.RTPS_parameter*self.sigma_b[species]) + ((1-self.RTPS_parameter)*sigma_a)
				ind0,ind1,ind2 = np.where(np.abs(sigma_a)>5e-16) #Machine precision can be a problem
				if len(ind0)>0:
					newoverold = sigma_RTPS[ind0,ind1,ind2]/sigma_a[ind0,ind1,ind2] #1D
					meanrebalance = np.mean(conc4D[:,ind0,ind1,ind2],axis=0)*(newoverold-1)
					conc4D[:,ind0,ind1,ind2] = (conc4D[:,ind0,ind1,ind2]*newoverold)-meanrebalance #Scale so sd is new_std and mean is old mean
					for i in self.ensemble_numbers: 
						self.gt[i].setSpecies3Dconc(species, conc4D[i-1,:,:,:])
	def saveRestartsAndScalingFactors(self,saveRestart=True, saveEmissions=True):
		for i in self.ensemble_numbers:
			if saveRestart:
				self.gt[i].saveRestart()
			if saveEmissions:
				self.gt[i].saveEmissions()
