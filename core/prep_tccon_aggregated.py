# aggregate TCCON observations and prepare hourly files for assimialtion in CHEEREIO (operator: tccon_tools) 
# prepared by Sina Voshtani 
# Modified by Drew Pendergrass (24 July 2025) for more flexible command line input 

#DEMO python prep_tccon_aggregated.py -i '/hpc/group/shindell/ap851/N2O/TCCON' -o '/hpc/group/shindell/ap851/N2O/TCCON4GC/2019' -s '2019-01-01T00:00:00' -e '2019-12-31T23:59:59' -s2k 'n2o'

import glob
import os
import numpy as np
import xarray as xr
import pandas as pd
import time
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser(description='Download GC input files.')
parser.add_argument('-i', '--input_path', type=str, help='Where are your input TCCON files, e.g. downloaded from tccondata.org')
parser.add_argument('-o', '--output_path', type=str, help='Where to save files.')
parser.add_argument('-p', '--input_file_pattern', type=str, default='*public.qc.nc', help='File pattern for input tccon files.')
parser.add_argument('-s', '--start_time', type=str, help='Start time, in format YYYY-MM-DDTHH:MM:SS. For example, 2023-01-01T00:00:00')
parser.add_argument('-e', '--end_time', type=str, help='End time, in format YYYY-MM-DDTHH:MM:SS. For example, 2023-01-01T00:00:00')
parser.add_argument('-s2k', '--species_to_keep', type=str, default='co,ch4,n2o,co2', help='Species to keep. If multiple, comma separated')
args = parser.parse_args()


# Specify the time range for TCCON aggregation 
start_time = np.datetime64(args.start_time)
end_time = np.datetime64(args.end_time)

# Specify the output directory
output_directory = args.output_path # path of outputs
file_pattern = args.input_file_pattern  
file_paths = glob.glob(os.path.join(args.input_path, file_pattern)) # provide the path of standard TCCON files (public and/or private files) - e.g., download from  https://tccondata.org/

# List of variables to keep
variables_to_keep = ['time', 'prior_time', 'prior_altitude', 'prior_pressure','prior_temperature', 'prior_h2o',  'ak_pressure',  'lat', 'long', 'pout']
# List of species-dependent variables and their corresponding error variables
species_dependent_vars = {}
# List of variables with two dimensions and their corresponding second dimensions
two_dimensional_vars = {}

species_to_keep = args.species_to_keep.split(',')
for species in species_to_keep:
    variables_to_keep.append(f'prior_{species}')
    variables_to_keep.append(f'ak_x{species}')
    variables_to_keep.append(f'x{species}')
    variables_to_keep.append(f'x{species}_error')
    species_dependent_vars[f'x{species}']=f'x{species}_error'
    species_dependent_vars[f'ak_x{species}']=f'x{species}_error'
    species_dependent_vars[f'prior_{species}']=f'x{species}_error'
    two_dimensional_vars[f'ak_x{species}']='ak_altitude'
    two_dimensional_vars[f'prior_{species}']='prior_altitude'


start_timestamp = time.time()

# Weighted median average function
def weighted_median(values, weights):
    sorted_indices = np.argsort(values)
    sorted_values = values[sorted_indices]
    sorted_weights = weights[sorted_indices]
    cumulative_weights = np.cumsum(sorted_weights)
    total_weight = np.sum(sorted_weights)
    median_index = np.searchsorted(cumulative_weights, total_weight / 2.0)
    weighted_median_value = sorted_values[median_index]
    return weighted_median_value

# Initialize an empty list to store datasets
datasets = []

for file_path in tqdm(file_paths, desc='Loading and Filtering Files', unit='file'):
    # Open the dataset
    ds = xr.open_dataset(file_path)

    # Extract the time variable
    time_variable = ds['time'].values.astype('datetime64[s]')

    # Check if the dataset has data within the specified time range
    if np.any((start_time <= time_variable) & (time_variable < end_time)):
        datasets.append(ds)

# Loop over each hour within the specified time range
total_hours = int((end_time - start_time) / np.timedelta64(1, 'h'))

for hour_offset in range(total_hours + 1):
    current_time = start_time + np.timedelta64(hour_offset, 'h')
    # Initialize an empty list to store datasets for the current hour
    hour_datasets = []

    # Loop over each loaded original dataset
    for ds in tqdm(datasets, desc=f'Processing Hour {current_time}', unit='file'):
        # Extract the time variable
        time_variable = ds['time'].values.astype('datetime64[s]')

        # Find indices corresponding to the current hour
        hour_indices = np.where((current_time <= time_variable) & (time_variable < current_time + np.timedelta64(1, 'h')))[0]

        # Continue to the next dataset if there is no data for the current hour
        if not hour_indices.any():
            continue

        # Create a subset with selected variables
        hour_subset_ds = ds[variables_to_keep].isel(time=hour_indices)
        
        #hour_subset_ds = hour_subset_ds.drop_duplicates('time', keep='first') #From Sina's code
        _, unique_index = np.unique(hour_subset_ds['time'], return_index=True) #Equivalent to Sina's code, but works with older versions of xr
        hour_subset_ds = hour_subset_ds.isel(time=unique_index)
        

        if hour_subset_ds.sizes['time'] == 0:
            print(f"No data within the specified time range for file: {file_path} and hour: {current_time}. Skipping to the next hour.")
        else:
            # Calculate hourly weighted median for species-dependent variables
            for species_var, error_var in species_dependent_vars.items():
                weights = 1 / (hour_subset_ds[error_var].values + 1e-10)  # Use respective error variable for weighting
                # Check if the variable has two dimensions
                if species_var in two_dimensional_vars:
                    second_dimension = two_dimensional_vars[species_var]
                    # Use isel to select relevant dimensions
                    values_2d = hour_subset_ds[species_var].values
                    weighted_median_result = np.apply_along_axis(weighted_median, axis=0,
                                                                 arr=values_2d,
                                                                 weights=weights)
                    # Ensure that the result has two dimensions
                    if len(weighted_median_result.shape) == 1:
                        # If 1D, reshape to 2D with the correct size
                        size_second_dim = len(hour_subset_ds[second_dimension])
                        weighted_median_result = weighted_median_result.reshape((1, size_second_dim))
                    # Create a DataArray with proper dimensions
                    hour_subset_ds[species_var] = xr.DataArray(
                        weighted_median_result,
                        dims=('time', second_dimension),
                        coords={'time': [hour_subset_ds['time'].values[0]], second_dimension: hour_subset_ds[second_dimension]}
                    )
                else:
                    hour_subset_ds[species_var] = np.apply_along_axis(weighted_median, axis=0,
                                                                     arr=hour_subset_ds[species_var].values,
                                                                     weights=weights)

            # Calculate hourly median for other variables
            hourly_median_ds = hour_subset_ds.resample(time='1H').median(dim='time')
            hourly_median_ds['time'] = hourly_median_ds['time'] + np.timedelta64(30, 'm')

            # Add a new variable time_utc
            hourly_median_ds['time_utc'] = hourly_median_ds['time'].dt.strftime('%Y-%m-%dT%H:%M:%S.%fZ')

            hourly_median_ds1_aligned = hourly_median_ds

            # Extract station location from the attribute 'short_location'
            station_location = ds.attrs.get('short_location', 'Unknown Location')

            # Use the station location as the station_id
            station_id_arr = xr.DataArray(
                np.array(station_location, dtype=str),
                dims=(),  # Empty tuple indicates a scalar dimension
            )

            hourly_median_ds1_aligned['station_id'] = station_id_arr

            # Append the aligned dataset to the list for the current hour
            hour_datasets.append(hourly_median_ds1_aligned)
            
    # Concatenate datasets along the latitude dimension only if there is data for the current hour
    if hour_datasets:
        combined_ds1 = xr.concat(hour_datasets, dim='site').squeeze('time')
        # Save the combined dataset to a new netCDF file for the current hour
        start_date_str = current_time.astype('M8[s]').astype('O').strftime('%Y%m%dT%H%M%S')
        end_date_str = (current_time + np.timedelta64(1, 'h')).astype('M8[s]').astype('O').strftime('%Y%m%dT%H%M%S')
        hour_combined_file_path = os.path.join(output_directory, f'tccon_avg_{start_date_str}_{end_date_str}.nc')
        combined_ds.to_netcdf(hour_combined_file_path)
        print(f'tccon hourly average saved to: {hour_combined_file_path}')

    # Move to the next hour
    current_time += np.timedelta64(1, 'h')

end_timestamp = time.time()
total_time = end_timestamp - start_timestamp
print(f'Total time taken: {total_time:.2f} seconds')

