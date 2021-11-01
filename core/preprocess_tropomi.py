#This file is modified from jacobian.py in the CH4_TROPOMI_INV repository, by M Sulprizio, 
#D Varon, L Shen and probably others not in the git history that I'm forgetting.
#See https://github.com/ACMG-CH4/CH4_TROPOMI_INV for the original

import glob
import numpy as np
import xarray as xr
import re
import pickle
import os
import pandas as pd 
import datetime
from shapely.geometry import Polygon

# Notes:
# ======
# - FIXED HARD-CODED PERTURBATION THING
# - The lat_ratio.csv file used for stratospheric correction is manually defined.
#   We may want to remove this feature entirely.
# - Not sure why we need both local times and UTC times. There are several
#   question-comments about this below.
# - Lu optimizes scaling factors, so he increases the sensitivities by a factor
#   of two, due to the 1.5 scaling perturbation. This is hard-coded and doesn't 
#   work if we are optimizing absolute emissions.
# - Not sure about some variables in the cal_weights() and remap()/remap2()
#   functions:
# 	    - data_type
# 	    - location
#	    - first_2
# - When Lu computes virtual TROPOMI column from GC data using the TROPOMI prior
#   and averaging kernel, he does it as the weighted mean mixing ratio [ppb] of
#   the relevant GC ground cells. Zhen does it as the weighted mean of number
#   of molecules instead. This requires saving out an additional GC diagnostic
#   variable -- something like the mass column in addition to PEDGE.
# - Need to double-check units of Jacobian [mixing ratio, unitless] vs units of
#   virtual TROPOMI column [ppb] in use_AK_to_GC().

# =============================================================================
#
#                                      Define functions
#
# =============================================================================

def save_obj(obj, name):
    """ Save something with Pickle. """
    
    with open(name , 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


# -----------------------------------------------------------------------------

def load_obj(name):
    """ Load something with Pickle. """

    with open( name, 'rb') as f:
        return pickle.load(f)


# -----------------------------------------------------------------------------
    
def read_tropomi(filename, species):
    """
    Read TROPOMI data and save important variables to dictionary.

    Arguments
        filename [str]  : TROPOMI netcdf data file to read
        species [str]   : Species string 

    Returns
        met      [dict] : Dictionary of important variables from TROPOMI:
	    	                - species
 			                - Latitude
			                - Longitude
       		                - QA value
 			                - UTC time
			                - Averaging kernel
			                - Prior profile
			                - Dry air subcolumns
			                - Latitude bounds
 			                - Longitude bounds
			                - Vertical pressure profile
    """

    # Initialize dictionary for TROPOMI data
    met = {}

    if filename=='TEST_NO2':
        filename = '/n/holylfs05/LABS/jacob_lab/dpendergrass/tropomi/raw/NO2/S5P_OFFL_L2__NO2____20211004T052433_20211004T070603_20598_02_020200_20211005T211455.nc'
    
    # Store species, QA, lat, lon, time, averaging kernel
    data = xr.open_dataset(filename, group='PRODUCT')
    if species=='NO2':
        met[species] = data['nitrogendioxide_tropospheric_column'].values[0,:,:]
    else:
        raise ValueError('Species not supported')

    met['qa_value'] = data['qa_value'].values[0,:,:]
    met['longitude'] = data['longitude'].values[0,:,:]
    met['latitude'] = data['latitude'].values[0,:,:]
    met['utctime'] = data['time_utc'].values[0,:]
    met['column_AK'] = data['averaging_kernel'].values[0,:,:,::-1]
    a = data['tm5_constant_a'].values[:,:]
    b = data['tm5_constant_b'].values[:,:]
    data.close()
    
    # Store methane prior profile, dry air subcolumns
    data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/INPUT_DATA')
    surface_pressure = data['surface_pressure'].values[0,:,:]/100                       # Pa -> hPa
    data.close()

    # Store lat, lon bounds for pixels
    data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/GEOLOCATIONS')
    met['longitude_bounds'] = data['longitude_bounds'].values[0,:,:,:]
    met['latitude_bounds'] = data['latitude_bounds'].values[0,:,:,:]
    data.close()

    #Pressure levels (hpa) with dimension level, latind, lonind,edge (bottom/top)
    pressures = np.zeros((np.shape(b)[0],np.shape(surface_pressure)[0],np.shape(surface_pressure)[1],np.shape(b)[1]))
    for i in range(np.shape(surface_pressure)[0]):
        for j in range(np.shape(surface_pressure)[1]):
            pressures[:,i,j,:]=a+(b*surface_pressure[i,j])

    met['pressures'] = pressures
    
    return met    

# -----------------------------------------------------------------------------

def cal_weights(Sat_p, GC_p):
    """
    Calculate pressure weights for TROPOMI & GEOS-Chem on a joint, uneven grid

    Arguments
        Sat_p   [float]    : Pressure edge from TROPOMI (13)          <--- 13-1 = 12 pressure levels
        GC_p    [float]    : Pressure edge from GEOS-Chem (48)        <--- 48-1 = 47 pressure levels

    Returns
        weights [dict]     : Pressure weights
    """

    # Combine Sat_p and GC_p into joint, uneven vertical grid
    Com_p = np.zeros(len(Sat_p)+len(GC_p))                     # <--- 61-1 = 60 pressure levels
    Com_p.fill(np.nan)
    data_type = np.zeros(len(Sat_p)+len(GC_p), dtype=int) 
    data_type.fill(-99)
    location = []
    i=0;j=0;k=0
    while ((i < len(Sat_p)) or (j < len(GC_p))):
        if i == len(Sat_p):
            Com_p[k] = GC_p[j]
            data_type[k] = 2
            j=j+1; k=k+1
            continue
        if j == len(GC_p):
            Com_p[k] = Sat_p[i]
            data_type[k] = 1  
            location.append(k)        
            i=i+1; k=k+1   
            continue
        if Sat_p[i] >= GC_p[j]:
            Com_p[k] = Sat_p[i]
            data_type[k] = 1        
            location.append(k)        
            i=i+1; k=k+1
        else:
            Com_p[k] = GC_p[j]
            data_type[k] = 2
            j=j+1; k=k+1
    
    # Find the first level of GC
    first_2 = -99
    for i in range(len(Sat_p)+len(GC_p)-1):
        if data_type[i] == 2:
            first_2 = i
            break    
    
    # Save data to dictionary
    weights = {}
    weights['data_type'] = data_type
    weights['Com_p'] = Com_p
    weights['location'] = location
    weights['first_2'] = first_2
    
    return weights


# -----------------------------------------------------------------------------

def remap(GC_CH4, data_type, Com_p, location, first_2):
    """
    Remap GEOS-Chem methane to TROPOMI vertical grid.

    Arguments
        GC_CH4    [float]   : Methane from GEOS-Chem (60) <---- 60 pressure levels?
        data_type [int]     : ****
        Com_p     [float]   : Combined TROPOMI + GEOS-Chem pressure levels
        location  [int]     : ****
        first_2   [int]     : ****

    Returns
        Sat_CH4   [float]   : GC methane in TROPOMI pressure coordinates
    """
    
    # ****
    conc = np.zeros(len(Com_p)-1,); conc.fill(np.nan)
    k=0
    for i in range(first_2, len(Com_p)-1):
        conc[i] = GC_CH4[k]
        if data_type[i+1] == 2:
            k=k+1
    if first_2 > 0:
        conc[:first_2] = conc[first_2]
    
    # Calculate the weighted mean methane for each layer
    delta_p = Com_p[:-1] - Com_p[1:]
    Sat_CH4 = np.zeros(12); Sat_CH4.fill(np.nan)
    for i in range(len(location)-1):
        start = location[i]
        end = location[i+1]
        fenzi = sum(conc[start:end]*delta_p[start:end])
        fenmu = sum(delta_p[start:end])
        Sat_CH4[i] = fenzi/fenmu   
    
    return Sat_CH4


# -----------------------------------------------------------------------------

def remap2(Sensi, data_type, Com_p, location, first_2):
    """
    Remap GEOS-Chem sensitivity data (from perturbation simulations) to TROPOMI vertical grid.

    Arguments
        Sensi     [float]   : 4D sensitivity data from GC perturbation runs, with dims (lon, lat, alt, cluster)****double-check dim order
        data_type [int]     : ****
        Com_p     [float]   : Combined TROPOMI + GEOS-Chem pressure levels
        location  [int]     : ****
        first_2   [int]     : ****

    Returns
        Sat_CH4   [float]   : GC methane in TROPOMI pressure coordinates
    """
    
    # ****
    MM = Sensi.shape[1]
    conc = np.zeros((len(Com_p)-1, MM))
    conc.fill(np.nan)
    k=0
    for i in range(first_2, len(Com_p)-1):
        conc[i,:] = Sensi[k,:]
        if data_type[i+1] == 2:
            k=k+1
    if first_2 > 0:
        conc[:first_2,:] = conc[first_2,:]
    
    # Calculate the weighted mean methane for each layer
    delta_p = Com_p[:-1] - Com_p[1:]
    delta_ps = np.transpose(np.tile(delta_p, (MM,1)))
    Sat_CH4 = np.zeros((12, MM)); Sat_CH4.fill(np.nan)
    for i in range(len(location)-1):
        start = location[i]
        end = location[i+1]
        fenzi = np.sum(conc[start:end,:]*delta_ps[start:end,:],0)
        fenmu = np.sum(delta_p[start:end])
        Sat_CH4[i,:] = fenzi/fenmu
    
    return Sat_CH4


# -----------------------------------------------------------------------------

def nearest_loc(loc_query, loc_grid, tolerance=0.5):
    """ Find the index of the nearest grid location to a query location, with some tolerance. """

    distances = np.abs(loc_grid - loc_query)
    ind = distances.argmin()
    if distances[ind] >= tolerance:
        return np.nan
    else:
        return ind


# -----------------------------------------------------------------------------


# =============================================================================
#
#                                      Run the code
#
# =============================================================================

if __name__ == '__main__':
    import sys

    startday = sys.argv[1]
    endday = sys.argv[2]
 
    # Reformat start and end days for datetime in configuration
    start = f'{startday[0:4]}-{startday[4:6]}-{startday[6:8]} 00:00:00'
    end = f'{endday[0:4]}-{endday[4:6]}-{endday[6:8]} 23:59:59'

    # Configuration
    correct_strato = False
    workdir = '.'
    Sensi_datadir = f'{workdir}/Sensi'
    Sat_datadir = f'{workdir}/data_TROPOMI'
    GC_datadir = f'{workdir}/data_GC'
    outputdir = f'{workdir}/data_converted'
    xlim = [-111,-95]
    ylim = [24,39]
    GC_startdate = np.datetime64(datetime.datetime.strptime(start, '%Y-%m-%d %H:%M:%S'))
    GC_enddate = np.datetime64(datetime.datetime.strptime(end, '%Y-%m-%d %H:%M:%S') - datetime.timedelta(days=1))
    print('Start:', start)
    print('End:', end)

    # Get TROPOMI data filenames for the desired date range
    print(Sat_datadir)
    allfiles = glob.glob(f'{Sat_datadir}/*.nc')
    Sat_files = []
    for index in range(len(allfiles)):
        filename = allfiles[index]
        shortname = re.split('\/', filename)[-1]
        shortname = re.split('\.', shortname)[0]
        strdate = re.split('\.|_+|T',shortname)[4]
        strdate2 = datetime.datetime.strptime(strdate, '%Y%m%d')
        if ((strdate2 >= GC_startdate) and (strdate2 <= GC_enddate)):
            Sat_files.append(filename)
    Sat_files.sort()
    print('Found', len(Sat_files), 'TROPOMI data files.')

    # Map GEOS-Chem to TROPOMI observation space
    # Also return Jacobian matrix if use_Sensi=True
    for filename in Sat_files:
        
        # Check if TROPOMI file has already been processed
        print('========================')
        temp = re.split('\/', filename)[-1]
        print(temp)
        date = re.split('\.',temp)[0]
        
        # If not yet processed, run use_AK_to_GC()
        if ~os.path.isfile(f'{outputdir}/{date}_GCtoTROPOMI.pkl'):
            if correct_strato:
                df = pd.read_csv('./lat_ratio.csv', index_col=0)
                lat_mid = df.index
                lat_ratio = df.values
                result = use_AK_to_GC(filename, n_clust, GC_startdate, GC_enddate, xlim, ylim, use_Sensi, Sensi_datadir, correct_strato, lat_mid, lat_ratio)
            else:
                print('Running use_AK_to_GC().')
                result = use_AK_to_GC(filename, n_clust, GC_startdate, GC_enddate, xlim, ylim, use_Sensi, Sensi_datadir)
        save_obj(result, f'{outputdir}/{date}_GCtoTROPOMI.pkl')

    print(f'Wrote files to {outputdir}')
