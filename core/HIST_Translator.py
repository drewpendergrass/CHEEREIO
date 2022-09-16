import xarray as xr
from glob import glob
import settings_interface as si 
from datetime import date,datetime,timedelta

#A class that takes history files and connects them with the main state vector and observation matrices
class HIST_Translator(object):
	def __init__(self, path_to_rundir,timeperiod,interval=None,verbose=1):
		self.verbose = verbose
		self.spc_config = si.getSpeciesConfig()
		self.hist_dir = f'{path_to_rundir}OutputDir'
		self.timeperiod = timeperiod
		self.interval = interval
	def globSubDir(self,timeperiod,useLevelEdge = False, useStateMet = False):
		subdir_lists = {}
		specconc_list = glob(f'{self.hist_dir}/GEOSChem.SpeciesConc*.nc4')
		specconc_list.sort()
		ts = [datetime.strptime(spc.split('.')[-2][0:13], "%Y%m%d_%H%M") for spc in specconc_list]
		if self.interval:
			specconc_list = [spc for spc,t in zip(specconc_list,ts) if (t>=timeperiod[0]) and (t<timeperiod[1]) and (t.hour % self.interval == 0)]
		else:
			specconc_list = [spc for spc,t in zip(specconc_list,ts) if (t>=timeperiod[0]) and (t<timeperiod[1])]
		subdir_lists['SpeciesConc'] = specconc_list
		if useStateMet:
			met_list = glob(f'{self.hist_dir}/GEOSChem.StateMet*.nc4')
			met_list.sort()
			met_ts = [datetime.strptime(met.split('.')[-2][0:13], "%Y%m%d_%H%M") for met in met_list]
			met_list = [met for met,t in zip(met_list,met_ts) if (t>=timeperiod[0]) and (t<timeperiod[1])]
			subdir_lists['StateMet'] = met_list
		if useLevelEdge:
			le_list = glob(f'{self.hist_dir}/GEOSChem.LevelEdgeDiags*.nc4')
			le_list.sort()
			le_ts = [datetime.strptime(le.split('.')[-2][0:13], "%Y%m%d_%H%M") for le in le_list]
			le_list = [le for le,t in zip(le_list,le_ts) if (t>=timeperiod[0]) and (t<timeperiod[1])]
			subdir_lists['LevelEdgeDiags'] = le_list
		return subdir_lists
	def combineHist(self,useLevelEdge=False, useStateMet = False):
		dataset=[]
		subdir_lists=self.globSubDir(self.timeperiod,useLevelEdge,useStateMet)
		for ind, specfile in subdir_lists['SpeciesConc']:
			hist_ds = xr.load_dataset(specfile)
			hist_vals = []
			le_vals = []
			met_vals = []
			for species in self.spc_config['HistorySpeciesConcToSave']:
				hist_vals.append(hist_ds[f'SpeciesConc_{species}'])
			if useLevelEdge:
				lefile = subdir_lists['LevelEdgeDiags'][ind]
				le_ds = xr.load_dataset(lefile)
				for lecoll in self.spc_config['HistoryLevelEdgeDiagsToSave']:
					le_vals.append(le_ds[lecoll])
			if useStateMet:
				metfile = subdir_lists['StateMet'][ind]
				met_ds = xr.load_dataset(metfile)
				for metcoll in self.spc_config['HistoryStateMetToSave']:
					met_vals.append(met_ds[metcoll])
			to_merge = hist_vals + le_vals + met_vals
			data_val = xr.merge(to_merge)
			dataset.append(data_val)
		dataset = xr.merge(dataset)
		return dataset
	def getArea(self):
		specconc_list=self.globSubDir(self.timeperiod,useLevelEdge=False,useStateMet=False)
		AREA = xr.load_dataset(specconc_list[0])[f'AREA']
		return AREA
