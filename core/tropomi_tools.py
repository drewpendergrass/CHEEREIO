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

def read_GC(date, use_Sensi=False, Sensi_datadir=None, correct_strato=False, lat_mid=None, lat_ratio=None):
    """
    Read GEOS-Chem data and save important variables to dictionary.

    Arguments
        date           [str]   : date of interest
        use_Sensi      [log]   : Are we trying to map GEOS-Chem sensitivities
                                 to TROPOMI observation space?
        Sensi_datadir  [str]   : If use_Sensi=True, this is the path to the GC
                                 sensitivity data
        correct_strato [log]   : Are we doing a latitudinal correction of
                                 GEOS-Chem stratospheric bias? 
        lat_mid        [float] : If correct_strato=True, this is the center
                                 latitude of each grid box
        lat_ratio      [float] : If correct_strato=True, this is the ratio for
                                 correcting GC stratospheric methane to match
                                 ACE-FTS
    
    Returns
        met            [dict]  : Dictionary of important variables from GEOS-Chem:
			                        - CH4
                                    - Latitude
                                    - Longitude
                                    - PEDGE
                                    - TROPP, if correct_strato=True
                                    - CH4_adjusted, if correct_strato=True
    """
    
    # Assemble file paths to GC output collections for input data
    month = int(date[4:6])
    file_species = f'GEOSChem.SpeciesConc.{date}00z.nc4'        
    file_pedge = f'GEOSChem.LevelEdgeDiags.{date}00z.nc4'    
    file_troppause = f'GEOSChem.StateMet.{date}00z.nc4'    
    
    # Read lat, lon, CH4 from the SpeciecConc collection
    filename = f'{GC_datadir}/{file_species}'
    data = xr.open_dataset(filename)
    LON = data['lon'].values
    LAT = data['lat'].values
    CH4 = data['SpeciesConc_CH4'].values[0,:,:,:]
    CH4 = CH4*1e9                                    # Convert to ppb
    CH4 = np.einsum('lij->jil', CH4)
    data.close()

    # Read PEDGE from the LevelEdgeDiags collection
    filename = f'{GC_datadir}/{file_pedge}'
    data = xr.open_dataset(filename)
    PEDGE = data['Met_PEDGE'].values[0,:,:,:]
    PEDGE = np.einsum('lij->jil', PEDGE)
    data.close()
    
    # If want to do latitudinal correction of stratospheric bias in GEOS-Chem:
    if correct_strato:
        
        # Read tropopause level from StateMet collection
        filename = f'{GC_datadir}/{file_troppause}'
        data = xr.open_dataset(filename)
        TROPP = data['Met_TropLev'].values[0,:,:]
        TROPP = np.einsum('ij->ji', TROPP)
        data.close()

        CH4_adjusted = CH4.copy()
        for i in range(len(LON)):
            for j in range(len(LAT)):
                l = int(TROPP[i,j])
                # Find the location of lat in lat_mid
                ind = np.where(lat_mid == LAT[j])[0][0]
                CH4_adjusted[i,j,l:] = CH4[i,j,l:]*lat_ratio[ind,month-1]
    
    # Store GC data in dictionary
    met = {}
    met['lon'] = LON
    met['lat'] = LAT
    met['PEDGE'] = PEDGE
    met['CH4'] = CH4
    if correct_strato:
        met['TROPP'] = TROPP
        met['CH4_adjusted'] = CH4_adjusted
    
    # If need to construct Jacobian, read sensitivity data from GEOS-Chem perturbation simulations
    if use_Sensi:
        filename = f'{Sensi_datadir}/Sensi_{date}.nc'
        data = xr.open_dataset(filename)
        Sensi = data['Sensitivities'].values
        Sensi = np.einsum('klji->ijlk', Sensi)
        data.close()
        # Now adjust the Sensitivity
        met['Sensi'] = Sensi

    return met


# -----------------------------------------------------------------------------

def read_all_GC(all_strdate, use_Sensi=False, Sensi_datadir=None, correct_strato=False, lat_mid=None, lat_ratio=None):
    """ 
    Call read_GC() for multiple dates in a loop. 

    Arguments
    Same as read_GC(), except instead of 'date' argument, use 
        allstr_date [list, str] : date strings 

    Returns
        met         [dict]      : Dictionary of dictionaries. Each
                                  sub-dictionary is returned by read_GC()
    """

    met={}
    for strdate in all_strdate:
        met[strdate] = read_GC(strdate, use_Sensi, Sensi_datadir, correct_strato, lat_mid, lat_ratio) 
    
    return met

# -----------------------------------------------------------------------------

def use_AK_to_GC(filename, n_clust, GC_startdate, GC_enddate, xlim, ylim, use_Sensi, Sensi_datadir, correct_strato=False, lat_mid=None, lat_ratio=None):
    """
    Map GEOS-Chem data to TROPOMI observation space.

    Arguments
        filename       [str]        : TROPOMI netcdf data file to read
        n_clust        [int]        : Number of clusters / state vector elements
        GC_startdate   [datetime64] : First day of inversion period, for GC and
                                      TROPOMI
        GC_enddate     [datetime64] : Last day of inversion period, for GC and
                                      TROPOMI
        xlim           [float]      : Longitude bounds for simulation domain
        ylim           [float]      : Latitude bounds for simulation domain
        use_Sensi      [log]        : Are we trying to map GEOS-Chem
                                      sensitivities to TROPOMI observation space?
        Sensi_datadir  [str]        : If use_Sensi=True, this is the path to the
                                      GC sensitivity data
        correct_strato [log]        : Are we doing a latitudinal correction of
                                      GEOS-Chem stratospheric bias?
        lat_mid        [float]      : If correct_strato=True, this is the center
                                      latitude of each grid box
        lat_ratio      [float]      : If correct_strato=True, this is the ratio
                                      for correcting GC stratospheric methane to
                                      match ACE-FTS

    Returns
        result         [dict]       : Dictionary with one or two fields:
 	   	 		                        - obs_GC : GEOS-Chem and TROPOMI methane data
	 					                            - TROPOMI methane
						                            - GC methane
						                            - TROPOMI lat, lon
						                            - TROPOMI lat index, lon index
				                      If use_Sensi=True, also include:
				                        - KK     : Jacobian matrix
    """
    
    # Read TROPOMI data
    TROPOMI = read_tropomi(filename)
    
    # We're only going to consider data within lat/lon/time bounds, and with QA > 0.5
    sat_ind = np.where((TROPOMI['longitude'] >  xlim[0])      & (TROPOMI['longitude'] <  xlim[1])     & 
                       (TROPOMI['latitude']  >  ylim[0])      & (TROPOMI['latitude']  <  ylim[1])     & 
                       (TROPOMI['localtime'] >= GC_startdate) & (TROPOMI['localtime'] <= GC_enddate)  &
                       (TROPOMI['qa_value']  >= 0.5))
    # [****Why are we using local times here? Shouldn't we be using UTC?]

    # Number of TROPOMI observations
    NN = len(sat_ind[0])
    print('Found', NN, 'TROPOMI observations.')

    # If need to build Jacobian from GC perturbation simulation sensitivity data:
    if use_Sensi:
        # Initialize Jacobian K
        temp_KK = np.zeros([NN,n_clust], dtype=np.float32)
        temp_KK.fill(np.nan)
    
    # Initialize array with NN rows, 6 columns: TROPOMI-CH4, GC-CH4, longitude, latitude, II, JJ
    temp_obs_GC = np.zeros([NN,6], dtype=np.float32)
    temp_obs_GC.fill(np.nan)
    
    # Initialize a list to store the dates we want to look at
    all_strdate = []

    # For each TROPOMI observation
    for iNN in range(NN):
        
        # Get the date and hour
        iSat = sat_ind[0][iNN] # lat index
        jSat = sat_ind[1][iNN] # lon index
        timeshift = int(TROPOMI['longitude'][iSat,jSat]/15*60)
        timeshift = 0 # Now I use UTC time instead of local time [****Why?]
        localtime = TROPOMI['utctime'][iSat] + np.timedelta64(timeshift, 'm') # local time
        localtime = pd.to_datetime(str(localtime))
        strdate = localtime.round('60min').strftime('%Y%m%d_%H')        
        all_strdate.append(strdate)
    all_strdate = list(set(all_strdate))

    # Read GEOS_Chem data for the dates of interest
    all_date_GC = read_all_GC(all_strdate, use_Sensi, Sensi_datadir, correct_strato, lat_mid, lat_ratio)
   
    # For each TROPOMI observation: 
    for iNN in range(NN):
        
        # Get GC data for the date of the observation:
        iSat = sat_ind[0][iNN]
        jSat = sat_ind[1][iNN]
        Sat_p = TROPOMI['pressures'][iSat,jSat,:]
        dry_air_subcolumns = TROPOMI['dry_air_subcolumns'][iSat,jSat,:]       # mol m-2
        priori = TROPOMI['methane_profile_apriori'][iSat,jSat,:]              # mol m-2
        AK = TROPOMI['column_AK'][iSat,jSat,:]
        timeshift = int(TROPOMI['longitude'][iSat,jSat]/15*60)
        timeshift = 0
        localtime = TROPOMI['utctime'][iSat]+np.timedelta64(timeshift,'m')    # local time
        localtime = pd.to_datetime(str(localtime))
        strdate = localtime.round('60min').strftime('%Y%m%d_%H')        
        GC = all_date_GC[strdate]
                
        # Find GC lats & lons closest to the corners of the TROPOMI pixel
        longitude_bounds = TROPOMI['longitude_bounds'][iSat,jSat,:]
        latitude_bounds = TROPOMI['latitude_bounds'][iSat,jSat,:]
        corners_lon = []
        corners_lat = []
        for k in range(4):
            iGC = nearest_loc(longitude_bounds[k], GC['lon'])
            jGC = nearest_loc(latitude_bounds[k], GC['lat'])
            corners_lon.append(iGC)
            corners_lat.append(jGC)
        GC_ij = [(x,y) for x in set(corners_lon) for y in set(corners_lat)]
        GC_grids = [(GC['lon'][i], GC['lat'][j]) for i,j in GC_ij]
        
        # Compute the overlapping area between the TROPOMI pixel and GC grid cells it touches
        overlap_area = np.zeros(len(GC_grids))
        dlon = GC['lon'][1]-GC['lon'][0]
        dlat = GC['lat'][1]-GC['lat'][0]
        # Polygon representing TROPOMI pixel
        p0 = Polygon(np.column_stack((longitude_bounds,latitude_bounds)))
        # For each GC grid cell that touches the TROPOMI pixel: 
        for ipixel in range(len(GC_grids)):
            # Define polygon representing the GC grid cell
            item = GC_grids[ipixel]
            ap1 = [item[0]-dlon/2, item[0]+dlon/2, item[0]+dlon/2, item[0]-dlon/2]
            ap2 = [item[1]-dlat/2, item[1]-dlat/2, item[1]+dlat/2, item[1]+dlat/2]        
            p2 = Polygon(np.column_stack((ap1, ap2)))
            # Calculate overlapping area as the intersection of the two polygons
            if p2.intersects(p0):
                  overlap_area[ipixel] = p0.intersection(p2).area
        
        # If there is no overlap between GC and TROPOMI, skip to next observation:
        if sum(overlap_area) == 0:
            continue                  

        # =======================================================
        #       Map GEOS-Chem to TROPOMI observation space
        # =======================================================
        
        # Otherwise, initialize tropomi virtual xch4 and virtual sensitivity as zero
        GC_base_posteri = 0   # virtual tropomi xch4
        GC_base_sensi = 0     # virtual tropomi sensitivity
        
        # For each GC ground cell that touches the TROPOMI pixel: 
        for ipixel in range(len(GC_grids)):
            
            # Get GC lat/lon indices for the cell
            iGC,jGC = GC_ij[ipixel]
            
            # Get GC pressure edges for the cell
            GC_p = GC['PEDGE'][iGC,jGC,:]
            
            # Get GC methane for the cell
            if correct_strato:
                GC_CH4 = GC['CH4_adjusted'][iGC,jGC,:]
            else:
                GC_CH4 = GC['CH4'][iGC,jGC,:]
                
            # Calculate pressure weights for the cell
            ww = cal_weights(Sat_p, GC_p)
            
            # Map GC methane to TROPOMI pressure levels
            Sat_CH4 = remap(GC_CH4, ww['data_type'], ww['Com_p'], ww['location'], ww['first_2'])   # ppb
            
            # Convert ppb to mol m-2
            Sat_CH4_2 = Sat_CH4 * 1e-9 * dry_air_subcolumns   # mol m-2
            
            # Derive the column-averaged XCH4 that TROPOMI would see over this ground cell
            tropomi_sees_ipixel = sum(priori + AK * (Sat_CH4_2-priori)) / sum(dry_air_subcolumns) * 1e9   # ppb
            
            # Weight by overlapping area (to be divided out later) and add to sum
            GC_base_posteri += overlap_area[ipixel] * tropomi_sees_ipixel  # ppb m2

            # If building Jacobian matrix from GC perturbation simulation
            # sensitivity data:
            if use_Sensi:
                
                # Get GC perturbation sensitivities at this lat/lon, for all vertical levels and clusters
                Sensi = GC['Sensi'][iGC,jGC,:,:]
                
                # Map the sensitivities to TROPOMI pressure levels
                Sat_CH4 = remap2(Sensi, ww['data_type'], ww['Com_p'], ww['location'], ww['first_2'])      # mixing ratio, unitless
                
                # Tile the TROPOMI averaging kernel
                AKs = np.transpose(np.tile(AK, (n_clust,1)))
                
                # Tile the TROPOMI dry air subcolumns
                dry_air_subcolumns_s = np.transpose(np.tile(dry_air_subcolumns, (n_clust,1)))   # mol m-2
                
                # Derive the change in column-averaged XCH4 that TROPOMI would see over this ground cell
                ap = np.sum(AKs*Sat_CH4*dry_air_subcolumns_s, 0) / sum(dry_air_subcolumns)   # mixing ratio, unitless
                
                # Weight by overlapping area (to be divided out later) and add to sum
                GC_base_sensi += overlap_area[ipixel] * ap  # m2

        # Compute virtual TROPOMI observation as weighted mean by overlapping area
        # i.e., need to divide out area [m2] from the previous step
        virtual_tropomi = GC_base_posteri/sum(overlap_area)
 
        # Save real and virtual TROPOMI data                        
        temp_obs_GC[iNN,0] = TROPOMI['methane'][iSat,jSat]     # TROPOMI methane
        temp_obs_GC[iNN,1] = virtual_tropomi                   # GC virtual TROPOMI methane
        temp_obs_GC[iNN,2] = TROPOMI['longitude'][iSat,jSat]   # TROPOMI longitude
        temp_obs_GC[iNN,3] = TROPOMI['latitude'][iSat,jSat]    # TROPOMI latitude 
        temp_obs_GC[iNN,4] = iSat                              # TROPOMI index of longitude
        temp_obs_GC[iNN,5] = jSat                              # TROPOMI index of lattitude

        if use_Sensi:
            # Compute TROPOMI sensitivity as weighted mean by overlapping area
            # i.e., need to divide out area [m2] from the previous step
            temp_KK[iNN,:] = GC_base_sensi/sum(overlap_area)
    
    # Output 
    result = {}

    # Always return the coincident TROPOMI and GC data
    result['obs_GC'] = temp_obs_GC

    # Optionally return the Jacobian
    if use_Sensi:
        result['KK'] = temp_KK
    
    return result
