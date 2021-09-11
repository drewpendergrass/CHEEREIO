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
	def __init__(self, nature_vals, nature_err_covariances,natureval_lats,natureval_lons,species_to_assimilate=None,testing=False):
		self.testing=testing
		#Create a dictionary for each species if multiple species
		if species_to_assimilate:
			self.values = dict(zip(species_to_assimilate, nature_vals))
			self.lats = dict(zip(species_to_assimilate,natureval_lats))
			self.lons = dict(zip(species_to_assimilate,natureval_lons))
			if type(nature_err_covariances) is np.ndarray:
				cov = []
				for i, nval in enumerate(nature_vals):
					cov.append(np.diag(nval*nature_err_covariances[i]))
					self.errs = dict(zip(species_to_assimilate, cov))
			else:
				self.errs = dict(zip(species_to_assimilate, species=None))
		else:
			self.values = nature_vals
			self.lats = natureval_lats
			self.lons = natureval_lons
			#If we're getting a matrix, that's the error covariances
			if type(nature_err_covariances) is np.ndarray:
				if np.shape(nature_err_covariances) != (len(nature_vals),len(nature_vals)):
					raise ValueError(f"Expected covariance matrix to be of dimension {len(nature_vals)}x{len(nature_vals)}; found dimension {np.shape(nature_err_covariances)}")
				else:
					self.errs = nature_err_covariances
			else:
				self.errs = np.diag(nature_vals*nature_err_covariances)
	def getObsVal(self,latind=None,lonind=None,species=None):
		if species:
			if latind:
				inds = self.getIndsOfInterest(latind,lonind,species)
				return self.values[species][inds]
			else:
				return self.values[species]
		else:
			if latind:
				inds = self.getIndsOfInterest(latind,lonind)
				return self.values[inds]
			else:
				return self.values
	def getObsErr(self,latind=None,lonind=None,species=None):
		if species:
			if latind:
				inds = self.getIndsOfInterest(latind,lonind,species)
				return self.errs[species][inds,inds]
			else:
				return self.errs[species]
		else:
			if latind:
				inds = self.getIndsOfInterest(latind,lonind)
				return self.errs[inds,inds]
			else:
				return self.errs
	def getIndsOfInterest(self,latind,lonind,species=None):
		data = tx.getSpeciesConfig(self.testing)
		loc_rad = float(data['LOCALIZATION_RADIUS_km'])
		gridlat,gridlon = tx.getLatLonVals(data,self.testing)
		latval = gridlat[latind]
		lonval = gridlon[lonind]
		if species:
			distvec = np.array([tx.calcDist_km(latval,lonval,a,b) for a,b in zip(self.lats[species],self.lons[species])])
		else:
			distvec = np.array([tx.calcDist_km(latval,lonval,a,b) for a,b in zip(self.lats,self.lons)])
		return np.where(distvec<=loc_rad)[0]
	def getObsLatLon(self,latind=None,lonind=None,species=None):
		latloninds = self.getIndsOfInterest(latind,lonind,species)
		if species:
			return [self.lats[species][latloninds],self.lons[species][latloninds]]
		else:
			return [self.lats[latloninds],self.lons[latloninds]]

#Parent class for all observation operators.
#Requires the full ensemble of concentrations, a 1D vector of "nature values",
#and a 2D vector of error covariance for those nature values.
#The key is the function H, which maps 3D concentrations in one ensemble member
#to a 1D vector of length nature_vals 
class ObsOperator(object):
	def __init__(self, nature_vals, nature_err_covariance,nature_lats,nature_lons,testing=False):
		self.obsinfo = ObservationInfo(nature_vals,nature_err_covariance,nature_lats,nature_lons,None,testing)
		self.testing=testing
	# Should map conc3D (all but last axis of conc4D) to 1D vector of shape nature_vals along with 1d vectors of lat and longitude indices to support localization.
	def H(self, conc3D,latinds=None,loninds=None):
		raise NotImplementedError
	def obsMeanAndPert(self, conc4D,latval=None,lonval=None):
		obsEns = np.zeros([len(self.obsinfo.getObsVal(latval,lonval)),np.shape(conc4D)[3]]) #Observation ensemble
		for i in range(np.shape(conc4D)[3]):
			if latval:
				latinds,loninds = tx.getIndsOfInterest(latval,lonval,self.testing)
				obsEns[:,i],_,_ = self.H(conc4D[:,:,:,i],latinds,loninds)
			else:
				obsEns[:,i],_,_ = self.H(conc4D[:,:,:,i])
		obsMean = np.expand_dims(np.mean(obsEns,axis = 1),axis=1)
		obsPert = np.transpose(np.transpose(obsEns)-np.transpose(obsMean))
		return [obsMean,obsPert]
	def obsDiff(self,vals,latval=None,lonval=None):
		return vals-self.obsinfo.getObsVal(latval,lonval)


#If we are doing an experiment where we model 'nature', this parent class wraps around
#the GC_Translator for the nature run. In particular, it will create an ObsOperator
#type class when needed. Needs to be given the error covariance matrix for the observations
#and the user needs to implement the H mapper.

class NatureHelper(object):
	def __init__(self, gt, species_to_assimilate, error_multipliers_or_matrices,testing=False):
		self.gt = gt
		self.testing = testing
		self.species_to_assimilate = species_to_assimilate
		nature_vecs = []
		nature_lats = []
		nature_lons = []
		for species in self.species_to_assimilate:
			conc3D = self.gt.getSpecies3Dconc(species)
			nature_vals,nature_lat,nature_lon = self.H(conc3D)
			nature_vecs.append(nature_vals)
			nature_lats.append(nature_lat)
			nature_lons.append(nature_lon)
		self.obs_info = ObservationInfo(nature_vecs,error_multipliers_or_matrices,nature_lats,nature_lons,self.species_to_assimilate,self.testing)
	def H(self, conc3D,latinds=None,loninds=None):
		raise NotImplementedError
	def getNatureVals(self,species,latind=None,lonind=None):
		return self.obs_info.getObsVal(latind,lonind,species)
	def getNatureErr(self,species,latind=None,lonind=None):
		return self.obs_info.getObsErr(latind,lonind,species)
	def getNatureLatLon(self,species,latind=None,lonind=None):
		return self.obs_info.getObsLatLon(latind,lonind,species)
	def makeObsOp(self,species,ObsOperatorClass,latind=None,lonind=None):
		nature_vals = self.getNatureVals(species,latind,lonind)
		nature_err_covariance = self.getNatureErr(species,latind,lonind)
		nature_lats,nature_lons = self.getObsLatLon(species,latind,lonind)
		return ObsOperatorClass(nature_vals,nature_err_covariance,nature_lats,nature_lons,None,self.testing)
	def makeR(self,latind=None,lonind=None):
		errmats = []
		for species in self.species_to_assimilate:
			errmats.append(self.getNatureErr(species,latind,lonind))
		R = block_diag(*errmats)
		return R


# ----------------------------------------------------------------- #
# -------------OPERATORS FOR TEST ASSIMILATION RUNS---------------- #
# ----------------------------------------------------------------- #

class SumOperator(ObsOperator):
	def H(self,conc3D,latinds=None,loninds=None):
		return column_sum(conc3D,latinds,loninds)

class SumNatureHelper(NatureHelper):
	def H(self,conc3D,latinds=None,loninds=None):
		return column_sum(conc3D,latinds,loninds)

# ----------------------------------------------------------------- #
# ------------------READY-MADE FUNCTIONS FOR H--------------------- #
# ----------------------------------------------------------------- #

#Sums up columns and optionally perturbs the result by a random normal variable 
#with mean given by bias and standard deviation given by err.
#If bias and/or error are None then the simple sum is returned.
#Returns a 1D array of the flattened 2D field.
#Takes numpy  array for species in question containing 3D concentrations.
def column_sum(DA_3d, latinds=None,loninds=None, bias=None, err=None):
	csum = np.sum(DA_3d,axis = 0)
	latvals,lonvals = tx.getLatLonVals()
	if latinds:
		csum = csum[latinds,loninds]
		latvals = latvals[latinds]
		lonvals = lonvals[loninds]
	if (bias is None) or (err is None):
		return [csum.flatten(),latvals,lonvals]
	else:
		csum += np.random.normal(bias, err, np.shape(csum))
		return [csum.flatten(),latvals,lonvals]

def surface_obs(DA_3d, latinds,loninds, bias=None, err=None):
	obs_vec = DA_3d[0,latinds,loninds]
	latvals,lonvals = tx.getLatLonVals()
	if (bias is None) or (err is None):
		return [obs_vec,latvals[latinds],lonvals[loninds]]
	else:
		obs_vec += np.random.normal(bias, err, np.shape(obs_vec))
		return [obs_vec,latvals[latinds],lonvals[loninds]]
