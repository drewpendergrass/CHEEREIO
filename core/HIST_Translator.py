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
		specconc_list = glob(f'{self.hist_dir}/GEOSChem.SpeciesConc*.nc4')
		specconc_list.sort()
		ts = [datetime.strptime(spc.split('.')[-2][0:13], "%Y%m%d_%H%M") for spc in specconc_list]
		if self.interval:
			specconc_list = [spc for spc,t in zip(specconc_list,ts) if (t>=timeperiod[0]) and (t<timeperiod[1]) and (t.hour % self.interval == 0)]
		else:
			specconc_list = [spc for spc,t in zip(specconc_list,ts) if (t>=timeperiod[0]) and (t<timeperiod[1])]
		if useStateMet:
			met_list = glob(f'{self.hist_dir}/GEOSChem.StateMet*.nc4')
			met_list.sort()
			met_ts = [datetime.strptime(met.split('.')[-2][0:13], "%Y%m%d_%H%M") for met in met_list]
			met_list = [met for met,t in zip(met_list,met_ts) if (t>=timeperiod[0]) and (t<timeperiod[1])]
		if useLevelEdge:
			le_list = glob(f'{self.hist_dir}/GEOSChem.LevelEdgeDiags*.nc4')
			le_list.sort()
			le_ts = [datetime.strptime(le.split('.')[-2][0:13], "%Y%m%d_%H%M") for le in le_list]
			le_list = [le for le,t in zip(le_list,le_ts) if (t>=timeperiod[0]) and (t<timeperiod[1])]
		if useStateMet and useLevelEdge:
			return [specconc_list,le_list,met_list]
		elif useStateMet:
			return [specconc_list,met_list]
		elif useLevelEdge:
			return [specconc_list,le_list]
		else:
			return specconc_list
	def combineHist(self,species,useLevelEdge=False, useStateMet = False):
		dataset=[]
		if useLevelEdge and useStateMet:
			specconc_list,le_list,met_list=self.globSubDir(self.timeperiod,useLevelEdge,useStateMet)
			for specfile,lefile,metfile in zip(specconc_list,le_list,met_list):
				hist_val = xr.load_dataset(specfile)[f'SpeciesConc_{species}']
				lev_val = xr.load_dataset(lefile)[f'Met_PEDGE']
				met_val = xr.load_dataset(metfile)[f'Met_AD']
				data_val = xr.merge([hist_val, lev_val,met_val])
				dataset.append(data_val)
		elif useLevelEdge:
			specconc_list,le_list=self.globSubDir(self.timeperiod,useLevelEdge,useStateMet)
			for specfile,lefile in zip(specconc_list,le_list):
				hist_val = xr.load_dataset(specfile)[f'SpeciesConc_{species}']
				lev_val = xr.load_dataset(lefile)[f'Met_PEDGE']
				data_val = xr.merge([hist_val, lev_val])
				dataset.append(data_val)
		elif useStateMet:
			specconc_list,met_list=self.globSubDir(self.timeperiod,useLevelEdge,useStateMet)
			for specfile,metfile in zip(specconc_list,met_list):
				hist_val = xr.load_dataset(specfile)[f'SpeciesConc_{species}']
				met_val = xr.load_dataset(metfile)[f'Met_AD']
				data_val = xr.merge([hist_val, met_val])
				dataset.append(data_val)
		else:
			specconc_list=self.globSubDir(self.timeperiod,useLevelEdge,useStateMet)
			for specfile in specconc_list:
				hist_val = xr.load_dataset(specfile)[f'SpeciesConc_{species}']
				dataset.append(hist_val)
		dataset = xr.merge(dataset)
		return dataset
	def getArea(self):
		specconc_list=self.globSubDir(self.timeperiod,useLevelEdge=False,useStateMet=False)
		AREA = xr.load_dataset(specconc_list[0])[f'AREA']
		return AREA
