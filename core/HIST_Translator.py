import xarray as xr
from glob import glob
import settings_interface as si 
from datetime import date,datetime,timedelta

#A class that takes history files and connects them with the main state vector and observation matrices
class HIST_Translator(object):
	def __init__(self, path_to_rundir,timeperiod,interval=None,verbose=1):
		self.verbose = verbose
		self.spc_config = si.getSpeciesConfig()
		gc_version = float(self.spc_config['GC_VERSION'][0:-2]) #major plus minor version
		if gc_version>=14.1:
			self.spcconc_name = "SpeciesConcVV"
		else:
			self.spcconc_name = "SpeciesConc" #Starting in 14.1 we have to specify VV
		self.hist_dir = f'{path_to_rundir}OutputDir'
		self.timeperiod = timeperiod
		self.interval = interval
	def makeCollDict(self,useLevelEdge = False, useStateMet = False, useObsPack = False):
		colls_to_grab = {'StateMet': {'use' : useStateMet, 'glob' : f'{self.hist_dir}/GEOSChem.StateMet*.nc4', 'diags' : 'HistoryStateMetToSave'},
		'LevelEdgeDiags': {'use' : useLevelEdge, 'glob' : f'{self.hist_dir}/GEOSChem.LevelEdgeDiags*.nc4', 'diags' : 'HistoryLevelEdgeDiagsToSave'},
		'ObsPack': {'use' : useObsPack, 'glob' : f'{self.hist_dir}/GEOSChem.ObsPack*.nc4','diags' : 'HistoryObsPackToSave'}
		}
		return colls_to_grab
	def globSubDir(self,timeperiod,useLevelEdge = False, useStateMet = False, useObsPack = False):
		subdir_lists = {}
		specconc_list = glob(f'{self.hist_dir}/GEOSChem.SpeciesConc*.nc4')
		specconc_list.sort()
		ts = [datetime.strptime(spc.split('.')[-2][0:13], "%Y%m%d_%H%M") for spc in specconc_list]
		if self.interval:
			specconc_list = [spc for spc,t in zip(specconc_list,ts) if (t>=timeperiod[0]) and (t<timeperiod[1]) and (t.hour % self.interval == 0)]
		else:
			specconc_list = [spc for spc,t in zip(specconc_list,ts) if (t>=timeperiod[0]) and (t<timeperiod[1])]
		subdir_lists['SpeciesConc'] = specconc_list
		#Grab additional collections which users might switch on
		colls_to_grab = self.makeCollDict(useLevelEdge, useStateMet, useObsPack)
		for coll in colls_to_grab:
			subdict = colls_to_grab[coll]
			if subdict['use']:
				met_list = glob(subdict['glob'])
				met_list.sort()
				met_ts = [datetime.strptime(met.split('.')[-2][0:13], "%Y%m%d_%H%M") for met in met_list]
				met_list = [met for met,t in zip(met_list,met_ts) if (t>=timeperiod[0]) and (t<timeperiod[1])]
				subdir_lists[coll] = met_list
		return subdir_lists
	def combineHist(self,useLevelEdge=False, useStateMet = False, useObsPack = False):
		dataset=[]
		to_merge=[]
		subdir_lists=self.globSubDir(self.timeperiod,useLevelEdge,useStateMet,useObsPack)
		colls_to_grab = self.makeCollDict(useLevelEdge, useStateMet, useObsPack)
		#Get species concentrations
		for ind, specfile in enumerate(subdir_lists['SpeciesConc']):
			hist_ds = xr.load_dataset(specfile)
			for species in self.spc_config['HistorySpeciesConcToSave']:
				to_merge.append(hist_ds[f'{self.spcconc_name}_{species}'])
		#Merge species concentrations and add to overall dataset
		data_val = xr.merge(to_merge)
		dataset.append(data_val)
		for coll in colls_to_grab:
			subdict = colls_to_grab[coll]
			if subdict['use']:
				#Obspack is a bit weirder and needs to be handled on its own
				if coll == "ObsPack":
					to_merge={}
					#Make dictionary entries, blank for now
					for lecoll in self.spc_config[subdict['diags']]:
						to_merge[lecoll] = []
					for ind in range(len(subdir_lists[coll])):
						lefile = subdir_lists[coll][ind]
						le_ds = xr.load_dataset(lefile)
						for lecoll in self.spc_config[subdict['diags']]:
							to_merge[lecoll].append(le_ds[lecoll])
					#Concatenate collections, merge doesn't like ragged arrays.
					for lecoll in self.spc_config[subdict['diags']]:
						if len(to_merge[lecoll]) > 0: #Not always obspack files
							data_val = xr.concat(to_merge[lecoll],'obs')
							dataset.append(data_val)
				else:
					to_merge=[]
					for ind in range(len(subdir_lists[coll])):
						lefile = subdir_lists[coll][ind]
						le_ds = xr.load_dataset(lefile)
						for lecoll in self.spc_config[subdict['diags']]:
							to_merge.append(le_ds[lecoll])
					data_val = xr.merge(to_merge)
					dataset.append(data_val)
		dataset = xr.merge(dataset)
		return dataset
	def reduceCombinedHistToSpecies(self,combinedHist,species):
		for spc in self.spc_config['HistorySpeciesConcToSave']:
			if spc != species:
				combinedHist = combinedHist.drop_vars(f'{self.spcconc_name}_{spc}')
		return combinedHist
	def getArea(self):
		specconc_list=self.globSubDir(self.timeperiod,useLevelEdge=False,useStateMet=False)
		AREA = xr.load_dataset(specconc_list[0])[f'AREA']
		return AREA
