import numpy as np
import xarray as xr
import toolbox as tx
from scipy.linalg import block_diag

#TODO: CHECK IF ERROR COVARIANCE CALCULATION IS REASONABLE

#Very simple class that wraps dictionaries linking (1) species name to 1D vector of observations, and
#(2) species name to 2D error covariance matrix.
#Basic constructor takes a list of species names, a list of 1D observation NumPy arrays, and 2D covariance NumPy arrays
#and pairs them.
#Class method also can produce diagonal covariance matrices if nature_err_covariances is a 1D numpy array containing percent errors
#Automatically subsets down according to lat and lon values.
class ObservationInfo(object):
	def __init__(self, species_to_assimilate, nature_vals, nature_err_covariances,natureval_lats,natureval_lons):
		self.values = dict(zip(species_to_assimilate, nature_vals))
		self.lats = dict(zip(species_to_assimilate,natureval_lats))
		self.lons = dict(zip(species_to_assimilate,natureval_lons))
		if type(nature_err_covariances) is np.ndarray:
			cov = []
			for i, nval in enumerate(nature_vals):
				cov.append(np.diag(nval*nature_err_covariances[i]))
				self.errs = dict(zip(species_to_assimilate, cov))
		else:
			self.errs = dict(zip(species_to_assimilate, nature_err_covariances))
	def getObsVal(self,species,latind=None,lonind=None):
		if latind:
			inds = self.getIndsOfInterest(species,latind,lonind)
			return self.values[species][inds]
		else:
			return self.values[species]
	def getObsErr(self,species,latind=None,lonind=None):
		if latind:
			inds = self.getIndsOfInterest(species,latind,lonind)
			return self.errs[species][inds,inds]
		else:
			return self.errs[species]
	def getIndsOfInterest(self,species,latind,lonind)
		data = tx.getSpeciesConfig()
		loc_rad = float(data['LOCALIZATION_RADIUS_km'])
		gridlat,gridlon = tx.getLatLonVals(data)
		latval = gridlat[latind]
		lonval = gridlon[lonind]
		distvec = np.array([tx.calcDist_km(latval,lonval,a,b) for a,b in zip(self.lats[species],self.lons[species])])
		return np.where(distvec<=loc_rad)[0]


#Parent class for all observation operators.
#Requires the full ensemble of concentrations, a 1D vector of "nature values",
#and a 2D vector of error covariance for those nature values.
#The key is the function H, which maps 3D concentrations in one ensemble member
#to a 1D vector of length nature_vals 
class ObsOperator(object):
	def __init__(self, nature_vals, nature_err_covariance):
		self.nature_vals = nature_vals
		self.nature_err_covariance = nature_err_covariance
	# Should map conc3D (all but last axis of conc4D) to 1D vector of shape nature_vals along with 1d vectors of lat and longitude indices to support localization.
	def H(self, conc3D):
		raise NotImplementedError
	def obsMeanAndPert(self, conc4D,latval=None,lonval=None):
		obsEns = np.zeros([len(self.nature_vals),np.shape(conc4D)[3]]) #Observation ensemble 
		for i in range(np.shape(conc4D)[3]):
			obsEns[:,i] = self.H(conc4D[:,:,:,i])
		obsMean = np.expand_dims(np.mean(obsEns,axis = 1),axis=1)
		obsPert = np.transpose(np.transpose(obsEns)-np.transpose(obsMean))
		return [obsMean,obsPert]
	def obsDiff(self,vals,latval=None,lonval=None):
		return vals-self.nature_vals


#If we are doing an experiment where we model 'nature', this parent class wraps around
#the GC_Translator for the nature run. In particular, it will create an ObsOperator
#type class when needed. Needs to be given the error covariance matrix for the observations
#and the user needs to implement the H mapper.

class NatureHelper(object):
	def __init__(self, gt, species_to_assimilate, error_multipliers_or_matrices):
		self.gt = gt
		self.species_to_assimilate = species_to_assimilate
		nature_vecs = []
		for species in self.species_to_assimilate:
			conc3D = self.gt.getSpecies3Dconc(species)
			nature_vals = self.H(conc3D)
			nature_vecs.append(nature_vals)
		self.obs_info = ObservationInfo(self.species_to_assimilate,nature_vecs,error_multipliers_or_matrices)
	def H(self, conc3D):
		raise NotImplementedError
	def getNatureVals(self,species):
		return self.obs_info.getObsVal(species)
	def getNatureErr(self,species):
		return self.obs_info.getObsErr(species)
	def makeObsOp(self,species,ObsOperatorClass,latind=None,lonind=None):
		nature_vals = self.getNatureVals(species)
		nature_err_covariance = self.getNatureErr(species)
		return ObsOperatorClass(nature_vals,nature_err_covariance)
	def makeR(self):
		errmats = []
		for species in self.species_to_assimilate:
			errmats.append(self.getNatureErr(species))
		R = block_diag(*errmats)
		return R


# ----------------------------------------------------------------- #
# -------------OPERATORS FOR TEST ASSIMILATION RUNS---------------- #
# ----------------------------------------------------------------- #

class SumOperator(ObsOperator):
	def H(self,conc3D):
		return column_sum(conc3D)

class SumNatureHelper(NatureHelper):
	def H(self,conc3D):
		return column_sum(conc3D)

# ----------------------------------------------------------------- #
# ------------------READY-MADE FUNCTIONS FOR H--------------------- #
# ----------------------------------------------------------------- #

#Sums up columns and optionally perturbs the result by a random normal variable 
#with mean given by bias and standard deviation given by err.
#If bias and/or error are None then the simple sum is returned.
#Returns a 1D array of the flattened 2D field.
#Takes numpy  array for species in question containing 3D concentrations.
def column_sum(DA_3d, bias=None, err=None):
	csum = np.sum(DA_3d,axis = 0)
	if (bias is None) or (err is None):
		return csum.flatten()
	else:
		csum += np.random.normal(bias, err, np.shape(csum))
		return csum.flatten()

def surface_obs(DA_3d, latinds,loninds, bias=None, err=None):
	obs_vec = DA_3d[0,latinds,loninds]
	if (bias is None) or (err is None):
		return obs_vec
	else:
		obs_vec += np.random.normal(bias, err, np.shape(obs_vec))
		return obs_vec
