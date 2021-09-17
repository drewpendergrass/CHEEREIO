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
		if self.testing:
			print("ObservationInfo constructor called.")
		#Create a dictionary for each species if multiple species
		if species_to_assimilate:
			#Error detection:
			if (len(nature_vals) != len(natureval_lats)) or (len(nature_vals) != len(natureval_lons)) or (len(natureval_lons) != len(natureval_lats)):
				raise ValueError('If passing lists of nature values/lats/lons, the lengths must match.')
			for i in range(len(nature_vals)):
				len_vals = len(nature_vals[i])
				len_lats = len(natureval_lats[i])
				len_lons = len(natureval_lons[i])
				if (len_vals != len_lats) or (len_vals != len_lons) or (len_lats != len_lons):
					raise ValueError(f'Lengths of vectors representing nature values ({len_vals}), latitudes ({len_lats}), and longitudes ({len_lons}) must match.') 
			self.values = dict(zip(species_to_assimilate.keys(), nature_vals))
			self.lats = dict(zip(species_to_assimilate.keys(),natureval_lats))
			self.lons = dict(zip(species_to_assimilate.keys(),natureval_lons))
			if type(nature_err_covariances) is np.ndarray:
				cov = []
				for i, nval in enumerate(nature_vals):
					cov.append(np.diag(nval*nature_err_covariances[i]))
					self.errs = dict(zip(species_to_assimilate.keys(), cov))
			else:
				raise NotImplementedError('Requires a list of error percentages')
		else:
			len_vals = len(nature_vals)
			len_lats = len(natureval_lats)
			len_lons = len(natureval_lons)
			if (len_vals != len_lats) or (len_vals != len_lons) or (len_lats != len_lons):
				raise ValueError(f'Lengths of vectors representing nature values ({len_vals}), latitudes ({len_lats}), and longitudes ({len_lons}) must match.') 
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
		if self.testing:
			print(f'ObservationInfo values have shape {np.shape(self.values)} and value {self.values}')
			print(f'ObservationInfo lats have shape {np.shape(self.lats)} and value {self.lats}')
			print(f'ObservationInfo lons have shape {np.shape(self.lons)} and value {self.lons}')
			print(f'ObservationInfo errs have shape {np.shape(self.errs)} and value {self.errs}')
			print(f'ObservationInfo constructor completed.')
	def getObsVal(self,latind=None,lonind=None,species=None):
		if self.testing:
			print(f'getObsVal called within ObservationInfo for species {species} and lat/lon inds {(latind,lonind)}.')
		if species:
			if not (latind is None):
				inds = self.getIndsOfInterest(latind,lonind,species)
				to_return = self.values[species][inds]
			else:
				to_return = self.values[species]
		else:
			if not (latind is None):
				inds = self.getIndsOfInterest(latind,lonind)
				to_return = self.values[inds]
			else:
				to_return = self.values
		if self.testing:
			print(f'getObsVal is returning array of shape {np.shape(to_return)}')
		return to_return
	def getObsErr(self,latind=None,lonind=None,species=None):
		if self.testing:
			print(f'getObsErr called within ObservationInfo for species {species} and lat/lon inds {(latind,lonind)}.')
		if species:
			if not (latind is None):
				inds = self.getIndsOfInterest(latind,lonind,species)
				to_return = np.zeros((len(inds),len(inds)))
				for i in range(len(inds)):
					to_return[i,:] = self.errs[species][inds[i],inds]
			else:
				to_return = self.errs[species]
		else:
			if not (latind is None):
				inds = self.getIndsOfInterest(latind,lonind)
				to_return = np.zeros((len(inds),len(inds)))
				for i in range(len(inds)):
					to_return[i,:] = self.errs[inds[i],inds]
			else:
				to_return = self.errs
		if self.testing:
			print(f'getObsErr is returning array of shape {np.shape(to_return)}')
		return to_return
	def getIndsOfInterest(self,latind,lonind,species=None):
		if self.testing:
			print(f'getIndsOfInterest called within ObservationInfo for species {species} and lat/lon inds {(latind,lonind)}.')
		data = tx.getSpeciesConfig(self.testing)
		loc_rad = float(data['LOCALIZATION_RADIUS_km'])
		origlat,origlon = tx.getLatLonVals(data,self.testing)
		latval = origlat[latind]
		lonval = origlon[lonind]
		if species:
			distvec = np.array([tx.calcDist_km(latval,lonval,a,b) for a,b in zip(self.lats[species],self.lons[species])])
		else:
			distvec = np.array([tx.calcDist_km(latval,lonval,a,b) for a,b in zip(self.lats,self.lons)])
		return np.where(distvec<=loc_rad)[0]
	def getObsLatLon(self,latind=None,lonind=None,species=None):
		if self.testing:
			print(f'getObsLatLon called within ObservationInfo for species {species} and lat/lon inds {(latind,lonind)}.')
		if species:
			if not (latind is None):
				latloninds = self.getIndsOfInterest(latind,lonind,species)
				to_return = [self.lats[species][latloninds],self.lons[species][latloninds]]
			else:
				to_return = [self.lats[species],self.lons[species]]
		else:
			if not (latind is None):
				latloninds = self.getIndsOfInterest(latind,lonind,species)
				to_return = [self.lats[latloninds],self.lons[latloninds]]
			else:
				to_return = [self.lats,self.lons]
		if self.testing:
			print(f'getObsLatLon is returning lat array of shape {np.shape(to_return[0])} and lon array of shape {np.shape(to_return[1])}')
		return to_return

#Parent class for all observation operators.
#Requires the full ensemble of concentrations, a 1D vector of "nature values",
#and a 2D vector of error covariance for those nature values.
#The key is the function H, which maps 3D concentrations in one ensemble member
#to a 1D vector of length nature_vals 
class ObsOperator(object):
	def __init__(self, nature_vals, nature_err_covariance,nature_lats,nature_lons,testing=False):
		self.testing=testing
		if self.testing:
			print('ObservationOperator constructor called.')
		self.obsinfo = ObservationInfo(nature_vals,nature_err_covariance,nature_lats,nature_lons,None,testing)
	# Should map conc3D (all but last axis of conc4D) to 1D vector of shape nature_vals along with 1d vectors of lat and longitude indices to support localization.
	def H(self, conc3D,latinds=None,loninds=None):
		raise NotImplementedError
	def obsMeanAndPert(self, conc4D,latval=None,lonval=None):
		if self.testing:
			print(f'obsMeanAndPert called within ObservationOperator for lat/lon inds {(latval,lonval)}.')
		obsEns = np.zeros([len(self.obsinfo.getObsVal(latval,lonval)),np.shape(conc4D)[3]]) #Observation ensemble
		for i in range(np.shape(conc4D)[3]):
			if not (latval is None):
				latinds,loninds = tx.getIndsOfInterest(latval,lonval,testing=self.testing)
				obsEns[:,i],_,_ = self.H(conc4D[:,:,:,i],latinds,loninds)
			else:
				obsEns[:,i],_,_ = self.H(conc4D[:,:,:,i],None,None)
		obsMean = np.mean(obsEns,axis = 1)
		obsPert = np.zeros(np.shape(obsEns))
		for j in range(np.shape(obsEns)):
			obsPert[:,j] = obsEns[:,j]-obsMean
		return [obsMean,obsPert]
	def obsDiff(self,vals,latval=None,lonval=None):
		return vals-self.obsinfo.getObsVal(latval,lonval)


#If we are doing an experiment where we model 'nature', this class wraps around
#the GC_Translator for the nature run. In particular, it will create an ObsOperator
#type class when needed. Needs to be given the error covariance matrix for the observations
#and species_to_assimilate is expected to be a dictionary mapping the user tag for this species to the species name 
class NatureHelper(object):
	def __init__(self, gt, species_to_assimilate, nature_h_functions, error_multipliers_or_matrices,testing=False):
		self.gt = gt
		self.testing = testing
		self.species_to_assimilate = species_to_assimilate
		nature_vecs = []
		nature_lats = []
		nature_lons = []
		if self.testing:
			print("NatureHelper constructor called.")
			data = tx.getSpeciesConfig(self.testing)
			obs_keys=list(data["OBSERVED_SPECIES"].keys())
			i=0
		for h,species in zip(nature_h_functions,self.species_to_assimilate.values()):
			conc3D = self.gt.getSpecies3Dconc(species)
			nature_vals,nature_lat,nature_lon = h(conc3D,testing=self.testing)
			nature_vecs.append(nature_vals)
			nature_lats.append(nature_lat)
			nature_lons.append(nature_lon)
			if self.testing:
				print(f"Applying function {data['NATURE_H_FUNCTIONS'][i]} to species {species} referenced by key {obs_keys[i]}.")
				i+=1
		self.obs_info = ObservationInfo(nature_vecs,error_multipliers_or_matrices,nature_lats,nature_lons,self.species_to_assimilate,self.testing)
		if self.testing:
			print("NatureHelper construction completed.")
	def getNatureVals(self,species,latind=None,lonind=None):
		if self.testing:
			print(f"Getting nature values for lat/lon inds {(latind,lonind)}.")
		return self.obs_info.getObsVal(latind,lonind,species)
	def getNatureErr(self,species,latind=None,lonind=None):
		if self.testing:
			print(f"Getting nature error for lat/lon inds {(latind,lonind)}.")
		return self.obs_info.getObsErr(latind,lonind,species)
	def getNatureLatLon(self,species,latind=None,lonind=None):
		if self.testing:
			print(f"Getting nature lat/lon for lat/lon inds {(latind,lonind)}.")
		return self.obs_info.getObsLatLon(latind,lonind,species)
	def makeObsOp(self,species,ObsOperatorClass):
		if self.testing:
			print(f"Making observation operator for species {species}.")
		nature_vals = self.getNatureVals(species)
		nature_err_covariance = self.getNatureErr(species)
		nature_lats,nature_lons = self.getNatureLatLon(species)
		return ObsOperatorClass(nature_vals,nature_err_covariance,nature_lats,nature_lons,self.testing)
	def makeR(self,latind=None,lonind=None):
		if self.testing:
			print(f"Making R for lat/lon inds {(latind,lonind)}.")
		errmats = []
		for species in self.species_to_assimilate:
			errmats.append(self.getNatureErr(species,latind,lonind))
		R = block_diag(*errmats)
		return R

def makeLatLonGrid(latvals,lonvals):
	latval_grid = np.transpose(np.tile(latvals,(len(lonvals),1)))
	lonval_grid = np.tile(lonvals,(len(latvals),1))
	return [latval_grid,lonval_grid]

# ----------------------------------------------------------------- #
# -------------OPERATORS FOR TEST ASSIMILATION RUNS---------------- #
# ----------------------------------------------------------------- #

class SumOperator(ObsOperator):
	def H(self,conc3D,latinds=None,loninds=None):
		return column_sum(conc3D,latinds,loninds,testing=self.testing)

class SurfaceOperator(ObsOperator):
	def H(self,conc3D,latinds,loninds):
		return surface_obs(conc3D,latinds,loninds,testing=self.testing)

# ----------------------------------------------------------------- #
# ------------------READY-MADE FUNCTIONS FOR H--------------------- #
# ----------------------------------------------------------------- #

#Sums up columns and optionally perturbs the result by a random normal variable 
#with mean given by bias and standard deviation given by err.
#If bias and/or error are None then the simple sum is returned.
#Returns a 1D array of the flattened 2D field.
#Takes numpy  array for species in question containing 3D concentrations.
def column_sum(DA_3d, latinds=None,loninds=None, bias=None, err=None, testing=False):
	if testing:
		print(f"Calculating column sum for lat/lon inds {(latinds,loninds)}.")
	csum = np.sum(DA_3d,axis = 0)
	latvals,lonvals = tx.getLatLonVals(testing=testing)
	latgrid,longrid = makeLatLonGrid(latvals,lonvals)
	if not (latinds is None):
		csum = csum[latinds,loninds]
		latvals = latgrid[latinds,loninds]
		lonvals = longrid[latinds,loninds]
	else:
		latvals = latgrid.flatten()
		lonvals = longrid.flatten()
	if (bias is None) or (err is None):
		return [csum.flatten(),latvals,lonvals]
	else:
		csum += np.random.normal(bias, err, np.shape(csum))
		return [csum.flatten(),latvals,lonvals]

def surface_obs(DA_3d, latinds=None,loninds=None, bias=None, err=None,testing=False):
	if testing:
		print(f"Calculating surface observations for lat/lon inds {(latinds,loninds)}.")
	latvals,lonvals = tx.getLatLonVals(testing=testing)
	latgrid,longrid = makeLatLonGrid(latvals,lonvals)
	if not (latinds is None):
		obs_vec = DA_3d[0,latinds,loninds]
		latvals = latgrid[latinds,loninds]
		lonvals = longrid[latinds,loninds]
	else:
		obs_vec = DA_3d[0,:,:].flatten()
		latvals = latgrid.flatten()
		lonvals = longrid.flatten()
	if (bias is None) or (err is None):
		return [obs_vec,latvals,lonvals]
	else:
		obs_vec += np.random.normal(bias, err, np.shape(obs_vec))
		return [obs_vec,latvals,lonvals]
