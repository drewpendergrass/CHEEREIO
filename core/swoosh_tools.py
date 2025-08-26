
from datetime import datetime
import xarray as xr
import numpy as np
import pandas as pd
import observation_operators as obsop
from scipy.interpolate import interp1d
import calendar

class SWOOSH_Translator(obsop.Observation_Translator):
	def __init__(self,verbose=1):
		super().__init__(verbose)
	def getObservations(self,specieskey,timeperiod, interval=None, includeObsError=False):
		species = self.spc_config['OBSERVED_SPECIES'][specieskey]
		swoosh = xr.open_dataset(self.spc_config["SWOOSH_dirs"][species]) #We assume the dir here is actually a file path
		if species == 'N2O':
			swoosh = swoosh['combinedn2oq'] #get N2O out
		else:
			raise ValueError(f'Species {species} not yet implemented for SWOOSH operator')
		#Do master filters the user expects for latitude
		if specieskey in list(self.spc_config["filter_obs_poleward_of_n_degrees"].keys()):
			filter_val = float(self.spc_config["filter_obs_poleward_of_n_degrees"][specieskey])
			if ~np.isnan(filter_val):
				swoosh = swoosh.where((swoosh.lat <= filter_val)&(swoosh.lat >= -1*filter_val), np.nan)
		#Filter by date
		idx = pd.DatetimeIndex(swoosh['time'].values)
		mon_start = idx.to_period('M').start_time
		mon_end   = idx.to_period('M').end_time
		start = pd.to_datetime(timeperiod[0])
		end   = pd.to_datetime(timeperiod[1])
		mask = (mon_end >= start) & (mon_start <= end)
		swoosh = swoosh.sel(time=mask)
		#Convert to PPB:
		swoosh*=1e3
		return swoosh
	def gcCompare(self,specieskey,swoosh,GC,GC_area=None,doErrCalc=True,useObserverError=False, prescribed_error=None,prescribed_error_type=None,transportError = None, errorCorr = None,minError=None):
		species = self.spc_config['OBSERVED_SPECIES'][specieskey]
		GC_SPC = obsop.getGC_SPC(GC,species,self.spc_config).resample(time="MS").mean()
		GC_P = GC['Met_PEDGE'].resample(time="MS").mean()
		latdiff = np.diff(swoosh.lat)/2
		latbins = np.concatenate([swoosh.lat[0:-1]-latdiff,[swoosh.lat[-1]-latdiff[-1],swoosh.lat[-1]+latdiff[-1]]])
		# Band means averaged over both lat & lon (per band)
		GC_SPC_latbands = GC_SPC.groupby_bins("lat", bins=latbins, right=False).mean(dim=("lat", "lon"))*1e9 #convert to ppb
		GC_P_latbands = GC_P.groupby_bins("lat", bins=latbins, right=False).mean(dim=("lat", "lon"))
		#Double check that we have equal dimensions in GC and SWOOSH. This can not be the case if we are on the edge of a month boundary
		if len(swoosh.time.values)!=len(GC_SPC_latbands.time.values):
			common = np.intersect1d(ym_code(GC_SPC_latbands), ym_code(swoosh))
			#Handle GC first
			mask = ym_code(GC_SPC_latbands).isin(common)
			GC_SPC_latbands = GC_SPC_latbands.where(mask, drop=True)
			GC_P_latbands = GC_P_latbands.where(mask, drop=True)
			#Handle swoosh now
			mask = ym_code(swoosh).isin(common)
			swoosh = swoosh.where(mask, drop=True)
		# Iterate through and interpolate GC
		gc_to_save = []
		obs_to_save = []
		lat_to_save = []
		lon_to_save = []
		levs_to_save = []
		time_to_save = []
		for t in range(len(swoosh.time.values)): #Sliced/averaged to have equal dimensions in GC and SWOOSH. This is monthly
			#Skip the month if less than days_allowed days present. Only need to check edge cases
			days_allowed=7
			if t==0:
				days_in_month = calendar.monthrange(GC.time[0].dt.year.values,GC.time[0].dt.month.values)[1]
				if days_in_month-GC.time[0].dt.day<=days_allowed: 
					continue #Skip, less than seven days present to end of month
			elif t==len(swoosh.time.values)-1:
				if GC.time[-1].dt.day<=days_allowed:
					continue #Skip, less than seven days present from start of month
			#If we have enough days, interpolate
			for i in range(len(swoosh.lat.values)): #Sliced/averaged to have equal dimensions in GC and SWOOSH
				gccol = GC_SPC_latbands.isel(time=t,lat_bins=i).values
				pcol = GC_P_latbands.isel(time=t,lat_bins=i).values
				pcol_mid = 0.5*(pcol[0:-1]+pcol[1:])
				f = interp1d(pcol_mid, gccol, kind='linear') #Interpolate GC to Swoosh levels
				GC_on_sat = f(swoosh.level.values)
				GC_lontile = np.tile(GC_on_sat, (len(GC.lon), 1)).flatten() #Tile (repeat) in lon dimension so it works with localization routines, then flatten to 1D
				gc_to_save.append(GC_lontile) #save flattened GC array
				obs = swoosh.isel(time=t,lat=i).values
				obs_lontile = np.tile(obs, (len(GC.lon), 1)).flatten()
				obs_to_save.append(obs_lontile)
				lat_to_save.append(np.repeat(swoosh.lat.values[i],len(GC_lontile))) #Save latitude. Everything is the same, so just repeat
				lon_to_save.append(np.repeat(GC.lon.values,len(swoosh.level.values))) #Repeats lon values each swoosh.levels in a row, as expected from flattening
				levs_to_save.append(np.tile(swoosh.level.values,len(GC.lon))) #Tiles lev values each lon times in a row, as expected from flattening
				time_to_save.append(np.repeat(swoosh.time.values[0],len(GC_lontile))) #Save time. Everything is the same, so just repeat
		#Do final quality checks
		df = pd.DataFrame()
		df['gc']=np.concatenate(gc_to_save)
		df['obs']=np.concatenate(obs_to_save)
		df['lats']=np.concatenate(lat_to_save)
		df['lons']=np.concatenate(lon_to_save)
		df['levs']=np.concatenate(levs_to_save)
		df['times']=np.concatenate(time_to_save)
		#Drop missing data
		df.dropna(inplace=True)
		#Drop data too low (i.e. high pressure) in MLS field
		df=df[df['levs']<=100]
		#Save out ObsData
		toreturn = obsop.ObsData(df['gc'].values,df['obs'].values,df['lats'].values,df['lons'].values,df['times'].values)
		toreturn.addData(level=df['levs'].values)
		return toreturn

def ym_code(da):
	return da.time.dt.year * 100 + da.time.dt.month

