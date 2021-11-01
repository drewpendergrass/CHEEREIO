## -------------------------------------------------------------------------##
## Load packages and set environment defaults
## -------------------------------------------------------------------------##
import numpy as np
import xarray as xr
import re
import pickle
import os
import sys
import pandas as pd
import datetime
import copy
import glob as glob

## -------------------------------------------------------------------------##
## Read in user preferences
## -------------------------------------------------------------------------##
# code_dir = '/n/home04/hnesser/TROPOMI_inversion/python'
# sat_data_dir = "/n/seasasfs02/CH4_inversion/InputData/Obs/TROPOMI/"
# GC_pressure_data_dir = "/n/holyscratch01/jacob_lab/hnesser/TROPOMI_inversion/jacobian_runs/TROPOMI_inversion_0000_final/OutputDir"
# GC_ch4_data_dir = "/n/holyscratch01/jacob_lab/hnesser/TROPOMI_inversion/jacobian_runs/TROPOMI_inversion_0034/OutputDir"
# output_dir = "/n/holyscratch01/jacob_lab/hnesser/TROPOMI_inversion/jacobian_runs/TROPOMI_inversion_0034/ProcessedDir/"
# MONTHS = [9, 10, 11, 12]
# jacobian = True

code_dir = sys.argv[1]
sat_data_dir = sys.argv[2]
GC_pressure_data_dir = sys.argv[3]
GC_ch4_data_dir = sys.argv[4]
output_dir = sys.argv[5]
MONTHS = [int(sys.argv[6])*4 -3 + i for i in range(4)]
jacobian = bool(sys.argv[7])
reprocess = False

# Load custom packages
sys.path.append(code_dir)
import gcpy as gc
import inversion_settings as s

print(f'Applying TROPOMI operator for {s.year}-{MONTHS}\n')

## -------------------------------------------------------------------------##
## Define functions
## -------------------------------------------------------------------------##
def save_obj(obj, name):
    with open(name , 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def nearest_loc(GC, TROPOMI):
    # Find the grid box and time indices corresponding to TROPOMI obs
    # i index
    iGC = np.abs(GC.lon.values.reshape((-1, 1))
                 - TROPOMI['longitude'].values.reshape((1, -1)))
    iGC = iGC.argmin(axis=0)

    # j index
    jGC = np.abs(GC.lat.values.reshape((-1, 1))
                 - TROPOMI['latitude'].values.reshape((1, -1)))
    jGC = jGC.argmin(axis=0)

    # Time index
    tGC = np.where(TROPOMI['utctime'].dt.hour == GC.time.dt.hour)[1]

    return iGC, jGC, tGC

def filter_tropomi(data, date, lon_min, lon_max, lat_min, lat_max):
    # Filter on qa_value
    data = data.where(data['qa_value'] > 0.5, drop=True)

    # Filter on masked methane
    data = data.where(data['xch4_corrected'] != 9.96921e36, drop=True)

    # Filter on lat/lon domain
    data = data.where((data['longitude_center'] >= lon_min) &
                      (data['longitude_center'] <= lon_max),
                      drop=True)
    data = data.where((data['latitude_center'] >= lat_min) &
                      (data['latitude_center'] <= lat_max),
                      drop=True)

    # Filter on dates
    data = data.where(((data['time'][:, 0] == int(date[:4]))
                       & (data['time'][:, 1] == int(date[4:6]))
                       & (data['time'][:, 2] == int(date[6:]))), drop=True)

    return data

def process_tropomi(data, date):
    # Do other processing
    # Add date variable
    dates = pd.DataFrame(data['time'].values[:, :-1],
                         columns=['year', 'month', 'day',
                                  'hour', 'minute', 'second'])
    dates = xr.DataArray(pd.to_datetime(dates),
                         dims=['nobs']).reset_index('nobs', drop=True)
    data = data.assign(utctime=dates)

    # Albedo and AOD have two columns [NIR, SWIR]. We select SWIR.
    # Both of these are needed for the albedo filter
    #data = data.where(data.nwin == 1, drop=True).squeeze()

    # Correct units from molecules/cm2 to mol/m2
    data['ch4_profile_apriori'] *= 1e4/6.02214e23
    data['dry_air_subcolumns'] *= 1e4/6.02214e23

    # Flip vertical direction
    data = data.sortby('nlayer', ascending=False)

    # Pressure information (hPa)
    pressure_interval = data['dp'].values.reshape((-1, 1))
    surface_pressure = data['surface_pressure'].values.reshape((-1, 1))

    # Create a pressures array corresponding to vertical levels
    # HN 2020/09/09 - converted from for loop to numpy
    z = data.nlayer.shape[0]
    pressures = (surface_pressure
                 - np.arange(z + 1).reshape((1, -1))*pressure_interval)
    pressures = xr.DataArray(pressures, dims=['nobs', 'ilayer'])
    pressures_mid = (surface_pressure
                     - (np.arange(z).reshape((1, -1)) + 0.5)*pressure_interval)
    pressures_mid = xr.DataArray(pressures_mid, dims=['nobs', 'nlayer'])
    data = data.assign(pressures=pressures, pressures_mid=pressures_mid)

    # Remove irrelevant variables
    data = data.drop(labels=['latitude_corners', 'longitude_corners',
                             'glintflag', 'altitude_levels',
                             'surface_altitude', 'time'])

    # Rename variables
    data = data.rename({'xch4_corrected' : 'methane',
                        'longitude_center' : 'longitude',
                        'latitude_center' : 'latitude',
                        'xch4_precision' : 'precision',
                        'surface_albedo' : 'albedo',
                        'aerosol_optical_thickness' : 'aerosol_optical_depth',
                        'xch4_column_averaging_kernel' : 'column_AK',
                        'ch4_profile_apriori' : 'methane_profile_apriori'})

    # Transpose
    data = data.transpose('nobs', 'nwin', 'nfov', 'nlayer', 'ilayer')

    return data

def get_diagnostic(data_dir, diag_name, date):
    short_date = date[:8]
    filename = os.path.join(data_dir,
                            'GEOSChem.'+diag_name+'.'+short_date+'_0000z.nc4')
    data = xr.open_dataset(filename)
    return data

def read_GC(ch4_data_dir, pedge_data_dir, date):
    # Start by downloading methane data (ppb)
    data = get_diagnostic(ch4_data_dir, 'SpeciesConc', date)
    data = data[['SpeciesConc_CH4']]*1e9

    # Get pressure information (hPa)
    pres = get_diagnostic(pedge_data_dir,
                          'LevelEdgeDiags', date)[['Met_PEDGE']]
    data = xr.merge([data, pres])
    pres.close()

    # Rename variables
    data = data.rename({'SpeciesConc_CH4' : 'CH4',
                        'Met_PEDGE' : 'PEDGE'})

    # Flip order
    data = data.transpose('time', 'lon', 'lat', 'lev', 'ilev')

    # Check that the data has all 24 hours
    if len(data.time) != 24:
        # print('GEOS-Chem Data does not contain 24 hours on %s.' % date)
        # print('Filling data with the first hour.')
        data = gc.fill_GC_first_hour(data)

    return data

def GC_to_sat_levels(GC_CH4, GC_edges, sat_edges):
    '''
    The provided edges for GEOS-Chem and the satellite should
    have dimension number of observations x number of edges
    '''
    # We want to account for the case when the GEOS-Chem surface
    # is above the satellite surface (altitude wise) or the GEOS-Chem
    # top is below the satellite top.. We do this by adjusting the
    # GEOS-Chem surface pressure up to the TROPOMI surface pressure
    idx_bottom = np.less(GC_edges[:, 0], sat_edges[:, 0])
    idx_top = np.greater(GC_edges[:, -1], sat_edges[:, -1])
    GC_edges[idx_bottom, 0] = sat_edges[idx_bottom, 0]
    GC_edges[idx_top, -1] = sat_edges[idx_top, -1]

    # Define vectors that give the "low" and "high" pressure
    # values for each GEOS-Chem and satellite layer.
    GC_lo = GC_edges[:, 1:][:, :, None]
    GC_hi = GC_edges[:, :-1][:, :, None]
    sat_lo = sat_edges[:, 1:][:, None, :]
    sat_hi = sat_edges[:, :-1][:, None, :]

    # Get the indices where the GC-to-satellite mapping, which is
    # a nobs x ngc x nsat matrix, is non-zero
    idx = (np.less_equal(sat_lo, GC_hi) & np.greater_equal(sat_hi, GC_lo))

    # Find the fraction of each GC level that contributes to each
    # TROPOMI level. We should first divide (to normalize) and then
    # multiply (to apply the map to the column) by the GC pressure
    # difference, but we exclude this (since it's the same as x1).
    GC_to_sat = np.minimum(sat_hi, GC_hi) - np.maximum(sat_lo, GC_lo)
    GC_to_sat[~idx] = 0

    # Now map the GC CH4 to the satellite levels
    GC_on_sat = (GC_to_sat*GC_CH4[:, :, None]).sum(axis=1)
    GC_on_sat = GC_on_sat/GC_to_sat.sum(axis=1)

    return GC_on_sat

def apply_avker(sat_avker, sat_prior, sat_pressure_weight, GC_CH4, filt=None):
    '''
    Apply the averaging kernel
    Inputs:
        sat_avker            The averaging kernel for the satellite
        sat_prior            The satellite prior profile in ppb
        sat_pressure_weight  The relative pressure weights for each level
        GC_CH4               The GC methane on the satellite levels
        filt                 A filter, optional
    '''
    if filt is None:
        filt = np.ones(sat_avker.shape[1])
    else:
        filt = filt.astype(int)

    GC_col = (filt*sat_pressure_weight
              *(sat_prior + sat_avker*(GC_CH4 - sat_prior)))
    GC_col = GC_col.sum(axis=1)
    return GC_col

## -------------------------------------------------------------------------##
## List all satellite files for the year and date defined
## -------------------------------------------------------------------------##
# List all raw netcdf TROPOMI files
raw_files = glob.glob(os.path.join(sat_data_dir, '*.nc'))
raw_files.sort()

# List all of the dates that have already been processed
if ~reprocess:
    saved_dates = glob.glob(os.path.join(output_dir, '*.pkl'))
    saved_dates = [d.split('/')[-1].split('_')[0] for d in saved_dates]
    saved_dates.sort()

# Create empty list
Sat_files = {}

# Iterate through the raw TROPOMI data
for index in range(len(raw_files)):
    filename = raw_files[index]

    # Get the date (YYYY, MM, and DD) of the raw TROPOMI file
    shortname = re.split('\/|\.', filename)[-2]
    strdate = re.split('_+|T', shortname)
    start_date = strdate[4]
    end_date = strdate[6]

    # Check if the start date is in range
    start = ((int(start_date[:4]) == s.year)
             and (int(start_date[4:6]) in MONTHS))

    # Check if the end date is in range
    end = ((int(end_date[:4]) == s.year)
           and (int(end_date[4:6]) in MONTHS))

    # Check if either date is in saved dates
    if ~reprocess:
        start = (start and (start_date not in saved_dates))
        end = (end and (end_date not in saved_dates))

    # Skip observations not in range
    if not (start or end):
        continue

    # Add the file to the list of Sat_files
    if start:
        if start_date in Sat_files.keys():
            Sat_files[start_date].append(filename)
        else:
            Sat_files[start_date] = [filename]
    elif end:
        if end_date in Sat_files.keys():
            Sat_files[end_date].append(filename)
        else:
            Sat_files[end_date] = [filename]

print('Number of dates: ', len(Sat_files))

## -------------------------------------------------------------------------##
## Apply the operator to each date of satellite observations
## -------------------------------------------------------------------------##
for date, filenames in Sat_files.items():
    # print('=========== %s ===========' % date)
    preprocess = lambda d: filter_tropomi(d, date,
                                          s.lon_min, s.lon_max,
                                          s.lat_min, s.lat_max)
    TROPOMI = xr.open_mfdataset(filenames, concat_dim='nobs',
                                combine='nested',
                                chunks=10000,
                                preprocess=preprocess)
    TROPOMI = process_tropomi(TROPOMI, date)

    if TROPOMI is None:
        print(f'{date} : 0 observations')
        # print('================================')
        continue

    # Get observation dimension (number of good observations in that single
    # observation file)
    NN = TROPOMI.nobs.shape[0]
    print(f'{date} : {NN} observations')

    # Then, read in the GC data for these dates.
    GC = read_GC(GC_ch4_data_dir, GC_pressure_data_dir, date)

    # Find the grid box and time indices corresponding to TROPOMI obs
    iGC, jGC, tGC = nearest_loc(GC, TROPOMI)

    # Then select GC accordingly
    GC_P = GC['PEDGE'].values[tGC, iGC, jGC, :]
    GC_CH4 = GC['CH4'].values[tGC, iGC, jGC, :]

    # And calculate TROPOMI XCH4 in ppb and the TROPOMI pressure weights
    TROP_CH4 = 1e9*(TROPOMI['methane_profile_apriori']
                    /TROPOMI['dry_air_subcolumns']).values
    TROP_PW = (-np.diff(TROPOMI['pressures'])/
               (TROPOMI['pressures'][:, 0] -
                TROPOMI['pressures'][:, -1]).values[:, None])

    # Map the GC CH4 onto the TROPOMI levels
    GC_on_sat = GC_to_sat_levels(GC_CH4, GC_P, TROPOMI['pressures'].values)

    # Finally, apply the averaging kernel
    GC_on_sat = apply_avker(TROPOMI['column_AK'].values,
                            TROP_CH4, TROP_PW, GC_on_sat)

    # Print out some general statistics....
    # print('Mean model - observation difference : ',
    #       np.mean(GC_on_sat - TROPOMI['methane'].values))
    # print('Standard deviation : ',
    #       np.std(GC_on_sat - TROPOMI['methane'].values))

    # Save out values
    # The columns are: OBS, MOD, LON, LAT, iGC, jGC, PRECISION,
    # ALBEDO_SWIR, ALBEDO_NIR, AOD, MOD_COL
    if jacobian:
        N = 10
    else:
        N = 2
    OBS_MOD = np.zeros([NN, N],dtype=np.float32)
    OBS_MOD[:, 0] = TROPOMI['methane']
    OBS_MOD[:, 1] = GC_on_sat
    if not jacobian:
        OBS_MOD[:, 2] = TROPOMI['longitude']
        OBS_MOD[:, 3] = TROPOMI['latitude']
        OBS_MOD[:, 4] = iGC
        OBS_MOD[:, 5] = jGC
        OBS_MOD[:, 6] = TROPOMI['precision']
        OBS_MOD[:, 7] = TROPOMI['albedo'][:,1]
        OBS_MOD[:, 8] = TROPOMI['albedo'][:,0]
        OBS_MOD[:, 9] = TROPOMI['aerosol_optical_depth'][:,1]

    save_obj(OBS_MOD, os.path.join(output_dir, date + '_GCtoTROPOMI.pkl'))
    # print('================================')

print('CODE FINISHED')
