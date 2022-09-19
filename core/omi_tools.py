from datetime import datetime,timedelta
from glob import glob
import pickle
import os.path
import xarray as xr
import numpy as np
import observation_operators as obsop
import functools

def read_omi(filename, species, filterinfo=None, includeObsError = False):
    """
    Read OMI data and save important variables to dictionary.

    Arguments
        filename [str]  : OMI netcdf data file to read

    Returns
        met   [dict] : Dictionary of important variables from OMI
    """

    if (filename=="test") and (species=="NO2"):
        filename = "/n/holylfs05/LABS/jacob_lab/dpendergrass/omi/NO2/2016/01/OMI-Aura_L2-OMNO2_2016m0116t2320-o61202_v003-2019m0819t154726.he5"

    # Initialize list for OMI data
    met = {}
    
    data = xr.open_dataset(filename, group='HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/')
    if species=="NO2":
        met['NO2'] = data['ColumnAmountNO2Trop'].values #Dimensions: time, Xtrack
        met['AmfTrop'] = data['AmfTrop'].values
        met['ScatteringWeight'] = data['ScatteringWeight'].values #Dimension: time, Xtrack, pressure level
        met['ScatteringWtPressure'] = data['ScatteringWtPressure'].values #Dimension: pressure level
        met['CloudRadianceFraction'] = data['CloudRadianceFraction'].values*0.001 # scale factor=0.001
        met['TerrainReflectivity'] = data['TerrainReflectivity'].values*0.001 # scale factor=0.001
        met['VcdQualityFlags'] = data['VcdQualityFlags'].values
        met['XTrackQualityFlags'] = data['XTrackQualityFlags'].values
        if includeObsError:
            met['Error'] = data['ColumnAmountNO2TropStd'].values #Dimensions: time, Xtrack
    data.close()

    data = xr.open_dataset(filename, group='HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/')

    met['longitude'] = data['Longitude'].values
    met['latitude'] = data['Latitude'].values

    utctime = np.zeros((len(data['Time'].values),1),dtype='object')
    tai93 = datetime(year=1993,month=1,day=1) # seconds since 1/1/1993 00:00 UTC
    for i in range(len(utctime)):
        utctime[i] = tai93 + timedelta(seconds=data['Time'].values[i])
    utctime = np.tile(utctime,(1,met['latitude'].shape[1]))
    met['utctime'] = utctime

    met['lat_corner_ll'] = data['FoV75CornerLatitude'].values[0,:,:]
    met['lat_corner_ul'] = data['FoV75CornerLatitude'].values[3,:,:]
    met['lat_corner_lr'] = data['FoV75CornerLatitude'].values[1,:,:]
    met['lat_corner_ur'] = data['FoV75CornerLatitude'].values[2,:,:]
    met['lon_corner_ll'] = data['FoV75CornerLongitude'].values[0,:,:]
    met['lon_corner_ul'] = data['FoV75CornerLongitude'].values[3,:,:]
    met['lon_corner_lr'] = data['FoV75CornerLongitude'].values[1,:,:]
    met['lon_corner_ur'] = data['FoV75CornerLongitude'].values[2,:,:]

    if species=="NO2":
        met['SolarZenithAngle'] = data['SolarZenithAngle'].values 

    data.close()

    met = clearEdgesFilterByQAAndFlatten(met)

    if filterinfo is not None:
        met = obsop.apply_filters(met,filterinfo)

    return met

#Clear swath edges, clear bad retrievals (bad QA), and flatten as CHEEREIO expects. Finally, drop nan rows. 
def clearEdgesFilterByQAAndFlatten(met):
    to_keep_by_flag = (met['VcdQualityFlags'] == 0) & ((met['XTrackQualityFlags'] == 0) | (met['XTrackQualityFlags'] == 255))
    met_toreturn = {}
    to_keep = []
    for key in met:
        #Skip flags; already incorporated
        if (key == "VcdQualityFlags") or (key == "XTrackQualityFlags"):
            continue
        temp = met[key]
        #Scattering weights have pressure dimension
        if key == 'ScatteringWeight':
            temp[:,0:5,:] = np.nan #remove swath edges
            temp[:,55:60,:] = np.nan
            met_toreturn[key] = temp[to_keep_by_flag]
        #Scattering weight pressure has only pressure dimension; just pass through
        elif key == 'ScatteringWtPressure':
            met_toreturn[key] = temp
        #Everything else handle as normal
        else:
            temp[:,0:5] = np.nan #remove swath edges
            temp[:,55:60] = np.nan
            met_toreturn[key] = temp[to_keep_by_flag]
            #Now we are going to drop Nans across the data; to this end we collect the places where there are no nans
            to_keep.append(~np.isnan(met_toreturn[key]))
    to_keep = functools.reduce(np.intersect1d, to_keep) #Where there are no nans across the data
    for key in met_toreturn:
        if key == "ScatteringWtPressure":
            continue
        if key == "ScatteringWeight":
            met_toreturn[key] = met_toreturn[key][to_keep,:] #Subset down to remove nans.
        else:
            met_toreturn[key] = met_toreturn[key][to_keep] #Subset down to remove nans.
    return met_toreturn
    

class OMI_Translator(obsop.Observation_Translator):
    def __init__(self,verbose=1):
        super().__init__(verbose)
    #Save dictionary of dates for later use
    def initialReadDate(self):
        sourcedirs = self.spc_config['OMI_dirs']
        OMI_date_dict = {}
        for key in list(sourcedirs.keys()):
            sourcedir = sourcedirs[key]
            obs_list = glob(f'{sourcedir}/**/*.he5', recursive=True)
            obs_list.sort()
            OMI_date_dict[key] = [datetime.strptime(obs[18:32], "%Ym%m%dt%H%M") for obs in obs_list]
        with open(f"{self.scratch}/omi_dates.pickle", 'wb') as handle:
            pickle.dump(OMI_date_dict, handle)
        return OMI_date_dict
    #Timeperiod is two datetime objects
    def globObs(self,species,timeperiod, interval=None):
        sourcedir = self.spc_config['OMI_dirs'][species]
        if os.path.exists(f"{self.scratch}/omi_dates.pickle"):
            with open(f"{self.scratch}/omi_dates.pickle", 'rb') as handle:
                OMI_date_dict = pickle.load(handle)
        else:
            OMI_date_dict = self.initialReadDate()
        obs_dates = OMI_date_dict[species]
        obs_list = glob(f'{sourcedir}/**/*.he5', recursive=True)
        obs_list.sort()
        if interval:
            obs_list = [obs for obs,t in zip(obs_list,obs_dates) if (t>=timeperiod[0]) and (t<timeperiod[1]) and ((t.hour % interval == 0) or (t.hour % interval == (interval-1)))]
        else:
            obs_list = [obs for obs,t in zip(obs_list,obs_dates) if (t>=timeperiod[0]) and (t<timeperiod[1])]
        return obs_list
    def getObservations(self,specieskey,timeperiod, interval=None, includeObsError=False):
        species = self.spc_config['OBSERVED_SPECIES'][specieskey]
        obs_list = self.globObs(species,timeperiod,interval)
        omi_obs = []
        filterinfo = {}
        if species=='NO2':
            if (self.spc_config['Extensions']['OMI_NO2']=="True") and (self.spc_config['OMI_NO2_FILTERS']=="True"): #Check first if extension is on before doing the OMI filtering
                filterinfo["OMI_NO2"] = [float(self.spc_config['OMI_NO2_filter_sza']),float(self.spc_config['OMI_NO2_filter_cloud_radiance_frac']),float(self.spc_config['OMI_NO2_filter_surface_albedo'])]
        if specieskey in list(self.spc_config["filter_obs_poleward_of_n_degrees"].keys()):
            filterinfo['MAIN']=[float(self.spc_config["filter_obs_poleward_of_n_degrees"][specieskey])]
        for obs in obs_list:
            omi_obs.append(read_omi(obs,species,filterinfo,includeObsError=includeObsError))
        met = {}
        for key in list(omi_obs[0].keys()):
            met[key] = np.concatenate([metval[key] for metval in omi_obs])
        return met
    def gcCompare(self,specieskey,TROPOMI,GC,GC_area=None,saveAlbedo=False,doErrCalc=True,useObserverError=False, prescribed_error=None,prescribed_error_type=None,transportError = None, errorCorr = None,minError=None):
        species = self.spc_config['OBSERVED_SPECIES'][specieskey]
        if species=='CH4':
            TROP_PRIOR = 1e9*(TROPOMI['methane_profile_apriori']/TROPOMI['dry_air_subcolumns'])
            synthetic_partial_columns = False
        elif species=='NO2':
            TROP_PRIOR=None
            synthetic_partial_columns = True
        TROP_PW = (-np.diff(TROPOMI['pressures'])/(TROPOMI['pressures'][:, 0] - TROPOMI['pressures'][:, -1])[:, None])
        returnStateMet = self.spc_config['SaveStateMet']=='True'
        if returnStateMet:
            GC_SPC,GC_P,GC_M,GC_area,i,j,t = obsop.getGCCols(GC,TROPOMI,species,returninds=True,returnStateMet=returnStateMet,GC_area=GC_area)
        else:
            GC_SPC,GC_P,GC_area,i,j,t = obsop.getGCCols(GC,TROPOMI,species,returninds=True,returnStateMet=returnStateMet,GC_area=GC_area)
        if species=='CH4':
            GC_SPC*=1e9 #scale to mol/mol
        #TROPOMI_ALL setting extension must be on!!
        memsetting = self.spc_config['LOW_MEMORY_TROPOMI_AVERAGING_KERNEL_CALC'] == 'True'
        if memsetting:
            batchsize = int(self.spc_config['LOW_MEMORY_TROPOMI_AVERAGING_KERNEL_BATCH_SIZE'])
        else:
            batchsize = None
        if returnStateMet:
            GC_on_sat,GC_M_on_sat = GC_to_sat_levels(GC_SPC, GC_P, TROPOMI['pressures'],GC_M = GC_M, lowmem=memsetting,batchsize=batchsize)
        else:
            GC_on_sat = GC_to_sat_levels(GC_SPC, GC_P, TROPOMI['pressures'],lowmem=memsetting,batchsize=batchsize)
            GC_M_on_sat = None
        GC_on_sat = apply_avker(TROPOMI['column_AK'],TROP_PW, GC_on_sat,TROP_PRIOR,GC_M_on_sat,GC_area)
        if self.spc_config['AV_TO_GC_GRID']=="True":
            superObsFunction = self.spc_config['SUPER_OBSERVATION_FUNCTION'][specieskey]
            additional_args_avgGC = {}
            if doErrCalc:
                if useObserverError:
                    additional_args_avgGC['obsInstrumentError'] = TROPOMI['Error']
                    additional_args_avgGC['modelTransportError'] = transportError
                elif prescribed_error is not None:
                    additional_args_avgGC['prescribed_error'] = prescribed_error
                    additional_args_avgGC['prescribed_error_type'] = prescribed_error_type
                if minError is not None:
                    additional_args_avgGC['minError'] = minError
                if errorCorr is not None:
                    additional_args_avgGC['errorCorr'] = errorCorr
            if saveAlbedo:
                additional_args_avgGC['albedo_swir'] = TROPOMI['albedo_swir']
                additional_args_avgGC['albedo_nir'] = TROPOMI['albedo_nir']
                additional_args_avgGC['blended_albedo'] = TROPOMI['blended_albedo']
            toreturn = obsop.averageByGC(i,j,t,GC,GC_on_sat,TROPOMI[species],doSuperObs=doErrCalc,superObsFunction=superObsFunction,**additional_args_avgGC)
        else:
            toreturn = obsop.ObsData(GC_on_sat,TROPOMI[species],TROPOMI['latitude'],TROPOMI['longitude'],TROPOMI['utctime'])
            if saveAlbedo:
                toreturn.addData(swir_av=TROPOMI['albedo_swir'],nir_av=TROPOMI['albedo_nir'],blended_av=TROPOMI['blended_albedo'])
            if doErrCalc and useObserverError:
                toreturn.addData(err_av=TROPOMI['Error'])
        return toreturn





