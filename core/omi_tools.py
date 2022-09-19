#THIS FILE IS OUT OF DATE AND BASED ON AN OLD VERSION OF CHEEREIO; MUST BE REWRITTEN BEFORE USE.

from datetime import datetime
from glob import glob
import pickle
import os.path
import xarray as xr
import numpy as np
import observation_operators as obsop

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
    
    data = xr.open_dataset(filename, group=group='HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/')
    if species=="NO2":
        met['NO2'] = data['ColumnAmountNO2Trop'].values #Dimensions: time, Xtrack
        met['AmfTrop'] = data['AmfTrop'].values
        met['ScatteringWeight'] = data['ScatteringWeight'].values
        met['ScatteringWtPressure'] = data['ScatteringWtPressure'].values
        met['CloudRadianceFraction'] = data['CloudRadianceFraction'].values*0.001 # scale factor=0.001
        met['TerrainReflectivity'] = data['TerrainReflectivity'].values*0.001 # scale factor=0.001
        met['VcdQualityFlags'] = data['VcdQualityFlags'].values
        if includeObsError:
            met['Error'] = data['ColumnAmountNO2TropStd'].values #Dimensions: time, Xtrack
    data.close()

    data = xr.open_dataset(filename, group=group='HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/')
    if species=="NO2":
        met['longitude'] = data['Longitude'].values
        met['latitude'] = data['Latitude'].values
        met['utctime'] = data['Time'].values 
        met['SolarZenithAngle'] = data['SolarZenithAngle'].values 
        met['FoV75CornerLongitude'] = data['FoV75CornerLongitude'].values 
        met['FoV75CornerLatitude'] = data['FoV75CornerLatitude'].values 
    data.close()

    if filterinfo is not None:
        met = obsop.apply_filters(met,filterinfo)

    return met


FROM VIRAL 

import numpy as np
import xarray as xr
import pandas as pd
import scipy.constants as const
import warnings
warnings.filterwarnings("ignore")

SEASdir='/n/home12/vshah/SEAS/'

# pressure to number density conversion factor
p_to_nd = ( 1.0 / const.R *
         10**2 * #hPa-> Pa
         const.value('Avogadro constant') * # mol/m3 -> molec/m3
         10**(-6)) # molec/m3 -> molec/cm3

# GEOS-Chem output at OMI crossing time (satellite diagnostic)
GC=xr.open_mfdataset(SEASdir+'BackgroundNO2/GEOSChem/rundirs/merra2_4x5_standard/'
                        'OutputDir_R20/ts_satellite.20150[678]*.bpch',engine='pseudonetcdf',
                          combine='by_coords')

# Seasonal mean
GC=GC.mean('time')

# Need to get tropopause height (could have archived it in the satellite diagnostic but will use the July mean here)
OutDir=SEASdir+'BackgroundNO2/GEOSChem/rundirs/merra2_4x5_standard/OutputDir_R0/'
Met=xr.open_dataset(OutDir+'GEOSChem.StateMet.20160701_0000z.nc4')

# Keep data only below the tropopause
tpp=Met.Met_TropP.isel(time=0).values
tpp3d=np.tile(tpp,(40,1,1))
pe=GC['PEDGE-$_PSURF'][1:]
pe=np.append(pe,pe[-2:-1,:,:],axis=0)
GC['IJ-AVG-$_NO2']=GC['IJ-AVG-$_NO2'].where(pe>tpp3d)

#select a subset of data & rename variables
GCdata=GC[['IJ-AVG-$_NO2','PEDGE-$_PSURF',"BXHGHT-$_BXHEIGHT",'DAO-3D-$_TMPU']]
GCdata=GCdata.rename({'IJ-AVG-$_NO2':'NO2MR','PEDGE-$_PSURF':'PE',
                     "BXHGHT-$_BXHEIGHT":'BXHGHT','DAO-3D-$_TMPU':'T'})

# GEOS-Chem pressure at mid-level
GC_pm=GCdata['PE']+np.diff(GCdata['PE'],axis=1,append=0)/2.0

GCdata['PM']=GCdata['PE']
GCdata.PM.values=GC_pm.values

# NO2 number density (molecules/cm3)
GCdata['NO2ND']=GCdata['NO2MR']*p_to_nd/GCdata['T']*GCdata['PM']*1e-9

# partial coumns (molecules/cm2)
GCdata['NO2vCol']=GCdata['NO2ND']*GCdata['BXHGHT']*1e2

# GC pressure mid points
press=GCdata.PM.mean(('latitude','longitude','time'))

# read scattering weights from binned OMI data and interpolate to GC pressure levels
# I had downloaded the OMI data for a season and binned it to the GEOS-Chem (4x5) grid
OMI=xr.open_dataset('./data/OMI2015_07-Copy1.nc').load()
sw=OMI.interp(press=press)["SW_bin"].mean('orbit').mean(('lat','lon'))

# GEOS-Chem VCD
GC_VCD=GCdata['NO2vCol'].mean(['latitude','longitude']).sum()
# GEOS-Chem SCD
GC_SCD=GCdata['NO2vCol'].mean(['latitude','longitude'])*sw).sum()

#AMF
GC_AMF = GC_SCD / GC_VCD

# You could either compare the GGEOS-Chem Slant Column Density (GC_SCD) with the OMI retrieved SCD,
# or GC_VCD with the OMI VCD re-calculated using the GEOS-Chem AMF
# OMI_VCD = OMI_SCD / GC_AMF





l2_dir = '/n/holylfs05/LABS/jacob_lab/nbalasus/omi/'
#west,east,south,north = 124,131,33,39
start_date = datetime(year=2021,month=4,day=21,hour=0,minute=0,second=0)
end_date = datetime(year=2021,month=4,day=21,hour=23,minute=59,second=0)
maxsza,maxcf = 70,0.3
data_fields = {'/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ColumnAmountNO2Trop':'column_amount',\
               '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ColumnAmountNO2TropStd':'column_uncertainty',\
               '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ScatteringWeight':'scattering_weights',\
               '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ScatteringWtPressure':'scattering_weights_pressure_level',\
               '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/CloudFraction':'cloud_fraction',\
               '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/VcdQualityFlags':'VCDQualityFlags',\
               '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/XTrackQualityFlags':'XTrackQualityFlags',\
               '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/FoV75CornerLatitude':'FoV75CornerLatitude',\
               '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/FoV75CornerLongitude':'FoV75CornerLongitude',\
               '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/SolarZenithAngle':'sza',\
               '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Latitude':'latc',\
               '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Longitude':'lonc',\
               '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Time':'time_utc'}

def F_subset_OMNO2(l2_dir,start_date,end_date,maxsza,maxcf,data_fields,west=None,east=None,south=None,north=None,verbose=True):

    # make list of all potential file names that could exist for our specified times
    l2_list = []

    # change wd to l2 files so I can use * for file path
    cwd = os.getcwd()
    os.chdir(l2_dir)

    days = (end_date-start_date).days+1
    dates = [start_date + timedelta(days=d) for d in range(days)]
    for date in dates:
        l2_list = l2_list + glob.glob('OMI-Aura_L2-OMNO2_'+date.strftime('%Ym%m%d')+'*.he5')

    # reset working directory to what is was before
    os.chdir(cwd)

    l2_data = {}

    for fn in l2_list:
        fn_dir = os.path.join(l2_dir,fn)
        if verbose == True:
            print('Loading '+fn)

        outp_he5 = {}
        with h5py.File(fn_dir,mode='r') as f:
            for product in data_fields:
                data = f[product]
                try:
                    scale_factor = data.attrs['ScaleFactor']
                    offset = data.attrs['Offset']
                    missing_value = data.attrs['MissingValue']

                except:
                    scale_factor = 1
                    offset = 0

                data = data[:]*scale_factor+offset

                outp_he5[data_fields[product]] = np.ma.masked_where(data == missing_value,data)

        # expand time to every pixel and convert to datetime
        time_utc = np.zeros((len(outp_he5['time_utc']),1),dtype='object')
        tai93 = datetime(year=1993,month=1,day=1) # seconds since 1/1/1993 00:00 UTC
        for i in range(len(time_utc)):
            time_utc[i] = tai93 + timedelta(seconds=outp_he5['time_utc'][i])
        outp_he5['time_utc'] = np.tile(time_utc,(1,outp_he5['latc'].shape[1]))

        # ll: lowerleft, ul: upperleft, ur: upperright, lr: lowerright
        outp_he5['lat_corner_ll'] = outp_he5['FoV75CornerLatitude'][0,:,:]
        outp_he5['lat_corner_ul'] = outp_he5['FoV75CornerLatitude'][3,:,:]
        outp_he5['lat_corner_lr'] = outp_he5['FoV75CornerLatitude'][1,:,:]
        outp_he5['lat_corner_ur'] = outp_he5['FoV75CornerLatitude'][2,:,:]
        outp_he5['lon_corner_ll'] = outp_he5['FoV75CornerLongitude'][3,:,:]
        outp_he5['lon_corner_ul'] = outp_he5['FoV75CornerLongitude'][0,:,:]
        outp_he5['lon_corner_lr'] = outp_he5['FoV75CornerLongitude'][1,:,:]
        outp_he5['lon_corner_ur'] = outp_he5['FoV75CornerLongitude'][2,:,:]

        f1 = outp_he5['sza'] <= maxsza
        f2 = outp_he5['cloud_fraction'] <= maxcf
        f3 = outp_he5['VCDQualityFlags'] == 0 & ((outp_he5['XTrackQualityFlags'] == 0) | (outp_he5['XTrackQualityFlags'] == 255))
        
        if north:
            f4 = np.min(np.array((outp_he5['lat_corner_ll'],outp_he5['lat_corner_ul'],outp_he5['lat_corner_ur'],outp_he5['lat_corner_lr'])),axis=0) >= south # south most lat
            f5 = np.max(np.array((outp_he5['lat_corner_ll'],outp_he5['lat_corner_ul'],outp_he5['lat_corner_ur'],outp_he5['lat_corner_lr'])),axis=0) <= north # north most lat
            f6 = np.min(np.array((outp_he5['lon_corner_ll'],outp_he5['lon_corner_ul'],outp_he5['lon_corner_ur'],outp_he5['lon_corner_lr'])),axis=0) >= west  # west most long
            f7 = np.max(np.array((outp_he5['lon_corner_ll'],outp_he5['lon_corner_ul'],outp_he5['lon_corner_ur'],outp_he5['lon_corner_lr'])),axis=0) <= east  # east most lon
        
        f8 = outp_he5['time_utc'] >= start_date
        f9 = outp_he5['time_utc'] <= end_date
        f10 = ~np.ma.getmask(outp_he5['column_amount']) & ~np.ma.getmask(outp_he5['column_uncertainty'])

        if north:
            validmask = f1 & f2 & f3 & f4 & f5 & f6 & f7 & f8 & f9 & f10
        else:
            validmask = f1 & f2 & f3 & f8 & f9 & f10

        if verbose == True:
            print('You have '+'%s'%np.sum(validmask)+' valid L2 pixels')

        # adding all of the other valid pixels to l2_data
        # don't want to have lat and lon bounds because I redefined them as corners
        for key in outp_he5.keys():
            if key not in ['FoV75CornerLatitude','FoV75CornerLongitude']:
                if key in l2_data:
                    l2_data[key] = np.append(l2_data[key],outp_he5[key][validmask])
                else:
                    l2_data[key] = outp_he5[key][validmask]
                    
    return l2_data

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
                filterinfo["OMI_NO2"] = [float(self.spc_config['TROPOMI_CH4_filter_blended_albedo']),float(self.spc_config['TROPOMI_CH4_filter_swir_albedo_low']),float(self.spc_config['TROPOMI_CH4_filter_swir_albedo_high']),float(self.spc_config['TROPOMI_CH4_filter_winter_lat']),float(self.spc_config['TROPOMI_CH4_filter_roughness']),float(self.spc_config['TROPOMI_CH4_filter_swir_aot'])]
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





