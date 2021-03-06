#Some of this code is based on work from Nick Balasus and Kang Sun based on https://amt.copernicus.org/articles/11/6679/2018/amt-11-6679-2018.pdf

import os
from netCDF4 import Dataset
from datetime import datetime,timedelta,time
import numpy as np
import glob
import h5py

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

class OMI_Translator(object):
    def __init__(self,testing=False):
        self.testing = testing
        self.spc_config = tx.getSpeciesConfig(self.testing)
        self.scratch = f"{self.spc_config['MY_PATH']}/{self.spc_config['RUN_NAME']}/scratch"
    #Save dictionary of dates for later use
    def initialReadDate(self):
        sourcedirs = self.spc_config['OMI_dirs']
        OMI_date_dict = {}
        for key in list(sourcedirs.keys()):
            sourcedir = sourcedirs[key]
            obs_list = glob(f'{sourcedir}/**/*.he5', recursive=True)
            obs_list.sort()2016m0809t0715
            OMI_date_dict[key] = [datetime.strptime(obs[18:32], "%Ym%m%dt%H%M") for obs in obs_list]
        with open(f"{self.scratch}/omi_dates.pickle", 'wb') as handle:
            pickle.dump(OMI_date_dict, handle)
        return OMI_date_dict
    #Timeperiod is two datetime objects
    def globObs(self,species,timeperiod):
        sourcedir = self.spc_config['OMI_dirs'][species]
        if os.path.exists(f"{self.scratch}/omi_dates.pickle"):
            with open(f"{self.scratch}/omi_dates.pickle", 'rb') as handle:
                OMI_date_dict = pickle.load(handle)
        else:
            OMI_date_dict = self.initialReadDate()
        obs_dates = OMI_date_dict[species]
        obs_list = glob(f'{sourcedir}/**/*.he5', recursive=True)
        obs_list.sort()
        obs_list = [obs for obs,t in zip(obs_list,obs_dates) if (t>=timeperiod[0]) and (t<timeperiod[1])]
        return obs_list
    def getSatellite(self,species,timeperiod):
        obs_list = self.globObs(species,timeperiod)
        trop_obs = []
        for obs in obs_list:
            trop_obs.append(read_tropomi(obs,species))
        met = {}
        for key in list(trop_obs[0].keys()):
            met[key] = np.concatenate([metval[key] for metval in trop_obs])
        return met
    def gcCompare(self,species,timeperiod,TROPOMI,GC):
        TROP_CH4 = 1e9*(TROPOMI['methane_profile_apriori']/TROPOMI['dry_air_subcolumns'])
        TROP_PW = (-np.diff(TROPOMI['pressures'])/(TROPOMI['pressures'][:, 0] - TROPOMI['pressures'][:, -1])[:, None])
        GC_CH4,GC_P = getGCCols(GC,TROPOMI,species)
        if species=='CH4':
            GC_CH4*=1e9 #scale to ppb
        GC_on_sat = GC_to_sat_levels(GC_CH4, GC_P, TROPOMI['pressures'])
        GC_on_sat = apply_avker(TROPOMI['column_AK'],TROP_CH4, TROP_PW, GC_on_sat)
        return [GC_on_sat,TROPOMI[species],TROPOMI['latitude'],TROPOMI['longitude']]



l2_data_omi = F_subset_OMNO2(l2_dir,start_date,end_date,maxsza,maxcf,data_fields,verbose=True)