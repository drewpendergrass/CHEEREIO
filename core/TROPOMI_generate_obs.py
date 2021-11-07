'''
This script generates the observation vector, the prior observation vector (i.e. F(xa)), and the observational error variances. It applies filters on albedo, latitude, and seaason to remove problematic TROPOMI observations. The variances are calcualted using the residual error method.

   **Inputs**

   | ----------------- | -------------------------------------------------- |
   | Input             | Description                                        |
   | ----------------- | -------------------------------------------------- |
   | prior_run         | A file or files containing the processed output of |
   |                   | the prior run of the forward model (after applying |
   |                   | the satellite averaging kernel. The input here can |
   |                   | be either a list of daily files or a single file   |
   |                   | containing all observations for the year.          |
   | ----------------- | -------------------------------------------------- |
   | filter_on_\       | A boolean flag indicating whether or not to filter |
   | blended_albedo    | on blended albedo, which should remove snow and    |
   |                   | ice covered scenes, as recommended by Lorente et   |
   |                   | al. 2021 and described in Wunch et al. 2011.       |
   | ----------------- | -------------------------------------------------- |
   | blended_albedo_ \ | The blended albedo threshold above which           |
   | threshold         | are removed from the observation vector. Lorente   |
   |                   | et al. find a value of 0.85 and Wunch et al. find  |
   |                   | a value of about 1. We use 1.1.                    |
   | ----------------- | -------------------------------------------------- |
   | filter_on_albedo  | A boolean flag indicating whether or not to filter |
   |                   | out scenes below the albedo_threshold, following   |
   |                   | the recommendation of de Gouw et al. 2020.         |
   | ----------------- | -------------------------------------------------- |
   | albedo_threshold  | The albedo threshold below which observations are  |
   |                   | removed. De Gouw et al. use 0.05. We do too.       |
   | ----------------- | -------------------------------------------------- |
   | filter_on_ \      | A boolean flag indicating whether or not to remove |
   | seaasonal_ \      | observations north of 50 degrees N during winter   |
   | latitude          | (DJF) months to further remove snow- and ice-      |
   |                   | covered scenes.                                    |
   | ----------------- | -------------------------------------------------- |
   | remove_ \         | A boolean flag indicating whether or not to remove |
   | latitudinal_bias  | the latitudinal bias in the model - observation    |
   |                   | difference with a first order polynomial.          |
   | ----------------- | -------------------------------------------------- |
   | analyze_biases    | A boolean flag indicating whether or not to        |
   |                   | analyze and plot spatial, latitudinal, seasonal,   |
   |                   | and albedinal, biases in the model - observation   |
   |                   | difference.                                        |
   | ----------------- | -------------------------------------------------- |
   | albedo_bins       | The albedo increments in which to bin the model -  |
   |                   | observation difference for statistical analysis.   |
   | ----------------- | -------------------------------------------------- |
   | lat_bins          | The latitude increments in which to bin the model  |
   |                   | - observation difference for statistical analysis. |
   | ----------------- | -------------------------------------------------- |
   | calculate_so      | A boolean flag indicating whether or not to        |
   |                   | calculate the observational error variances using  |
   |                   | the residual error method (Heald et al. 2004).     |
   | ----------------- | -------------------------------------------------- |

   **Outputs**

   | ----------------- | -------------------------------------------------- |
   | Output            | Description                                        |
   | ----------------- | -------------------------------------------------- |
   | y.nc              | The observation vector containing bias-corrected   |
   |                   | TROPOMI observations.                              |
   | ----------------- | -------------------------------------------------- |
   | ya.nc             | The prior observation vector containing the output |
   |                   | of the prior simulation, i.e. F(xa).               |
   | ----------------- | -------------------------------------------------- |
   | so.nc             | The observational error variances calculated using |
   |                   | the residual error method.                         |
   | ----------------- | -------------------------------------------------- |
   | blended_albedo_ \ | Monthly plots of observed XCH4 vs. blended albedo. |
   | filter_mm*.png    |                                                    |
   | ----------------- | -------------------------------------------------- |
   | prior_bias*.png   | Scatter plot of model vs. observation.             |
   | ----------------- | -------------------------------------------------- |
   | prior_spatial_ \  | Map of the model - observation bias on the inverse |
   | bias*.png         | grid.                                              |
   | ----------------- | -------------------------------------------------- |
   | prior_ \          | Plot of the mean and standard deviation of the     |
   | latitudinal_ \    | model - observation bias binned by lat_bin.        |
   | bias*.png         |                                                    |
   | ----------------- | -------------------------------------------------- |
   | prior_seasonal \  | Plot of the mean and standard deviation of the     |
   | bias*.png         | model - observation bias binned by month.          |
   | ----------------- | -------------------------------------------------- |
   | prior_seasonal \  | Plot of the mean and standard deviation of the     |
   | latitudinal_ \    | model - observation bias binned by lat_bin and     |
   | bias*.png         | month.                                             |
   | ----------------- | -------------------------------------------------- |
   | prior_albedo_ \   | Plot of the mean and standard deviation of the     |
   | bias*.png         | model - observation bias binned by albedo_bin.     |
   | ----------------- | -------------------------------------------------- |
   | observations*.png | Seasonal plots of the observations averaged over   |
   |                   | the inversion grid.                                |
   | ----------------- | -------------------------------------------------- |
   | errors*.png       | Seaasonal plots of the standard deviation averaged |
   |                   | over the inversion grid.                           |
   | ----------------- | -------------------------------------------------- |
'''

from os.path import join
from os import listdir
import sys
import copy
import calendar as cal
import xarray as xr
import numpy as np
from numpy.polynomial import polynomial as p
import pandas as pd
import matplotlib.pyplot as plt
pd.set_option('display.max_columns', 10)

## ------------------------------------------------------------------------ ##
## Set user preferences
## ------------------------------------------------------------------------ ##
# # Local preferences
# base_dir = '/Users/hannahnesser/Documents/Harvard/Research/TROPOMI_Inversion/'
# code_dir = base_dir + 'python'
# data_dir = base_dir + 'inversion_data'
# output_dir = base_dir + 'inversion_data'
# plot_dir = base_dir + 'plots'

# Cannon long-path preferences
base_dir = '/n/holyscratch01/jacob_lab/hnesser/TROPOMI_inversion/jacobian_runs/TROPOMI_inversion_0000/'
code_dir = '/n/home04/hnesser/TROPOMI_inversion/python'
data_dir = f'{base_dir}ProcessedDir'
output_dir = '/n/seasasfs02/hnesser/TROPOMI_inversion/inversion_data'

# Cannon preferences
# code_dir = sys.argv[1]
# base_dir = sys.argv[2]
# data_dir = f'{base_dir}ProcessedDir'
# output_dir = sys.argv[3]
# plot_dir = None

# Import Custom packages
sys.path.append(code_dir)
import config
config.SCALE = config.PRES_SCALE
config.BASE_WIDTH = config.PRES_WIDTH
config.BASE_HEIGHT = config.PRES_HEIGHT
import gcpy as gc
import troppy as tp
import format_plots as fp
import inversion_settings as settings

# The prior_run can either be a list of files or a single file
# with all of the data for simulation
prior_run = f'{settings.year}.pkl'
# prior_run = [f'{settings.year}{mm:02d}{dd:02d}_GCtoTROPOMI.pkl'
#              for mm in settings.months
#              for dd in settings.days]
# prior_run.sort()

# Define the blended albedo threshold
filter_on_blended_albedo = True
blended_albedo_threshold = 0.75

# Define a plain old albedo threshold
filter_on_albedo = True
albedo_threshold = 0.05

# Define a seasonal latitudinal filter
filter_on_seasonal_latitude = True

# Remove latitudinal bias
remove_latitudinal_bias = True
lat_bins = np.arange(10, 65, 5)

# Which analyses do you wish to perform?
analyze_biases = True
albedo_bins = np.arange(0, 1.1, 0.05)
lat_bins = np.arange(10, 65, 5)

# Calculate the error variances?
calculate_so = True

## ------------------------------------------------------------------------ ##
## Define functions
## ------------------------------------------------------------------------ ##
# Define an empty suffix for file names
suffix = ''
if filter_on_blended_albedo:
    suffix += '_BAF'
if filter_on_albedo:
    suffix += '_AF'
if filter_on_seasonal_latitude:
    suffix += '_LF'
if remove_latitudinal_bias:
    suffix += '_BC'

def apply_filter(data, criteria, filter_name):
    old_nobs = data.shape[0]
    data = data[criteria]
    new_nobs = data.shape[0]

    # Print information
    print(f'\nData is filtered on {filter_name}.')
    print(f'    {100*new_nobs/old_nobs:.1f}% of data is preserved.')
    print(f'    {new_nobs} observations remain.')
    # print(f'    Statistics:')
    # summ = data.groupby('MONTH')[['MOD', 'OBS', 'DIFF']]
    # summ = pd.concat([summ.min(), summ.max()], axis=1)
    # print(summ)

    return data

## ------------------------------------------------------------------------ ##
## Calculate y
## ------------------------------------------------------------------------ ##
if type(prior_run) == list:
    ## ----------------------------------------- ##
    ## Load data for the year
    ## ----------------------------------------- ##
    data = np.array([]).reshape(0, 12)
    for file in prior_run:
        # Check if that file is in the data directory
        if file not in listdir(data_dir):
            print(f'{file} is not in the data directory.')
            continue

        # Get the month
        month = int(file[4:6])
        day = int(file[6:8])

        # Load the data. The columns are: 0 OBS, 1 MOD, 2 LON, 3 LAT,
        # 4 iGC, 5 jGC, 6 PRECISION, 7 ALBEDO_SWIR, 8 ALBEDO_NIR, 9 AOD,
        # (15 10 total columns)
        new_data = gc.load_obj(join(data_dir, file))
        new_data = new_data[:, :10]
        new_data = np.insert(new_data, 10, month, axis=1) # add month
        new_data = np.insert(new_data, 11, day, axis=1)

        data = np.concatenate((data, new_data))

    ## ----------------------------------------- ##
    ## Basic data formatting
    ## ----------------------------------------- ##
    # Create a dataframe from the data
    columns = ['OBS', 'MOD', 'LON', 'LAT', 'iGC', 'jGC', 'PREC',
               'ALBEDO_SWIR', 'ALBEDO_NIR', 'AOD', 'MONTH', 'DAY']
    data = pd.DataFrame(data, columns=columns)

    # Calculate blended albedo
    data['BLENDED_ALBEDO'] = tp.blended_albedo(data,
                                               data['ALBEDO_SWIR'],
                                               data['ALBEDO_NIR'])

    # Subset data
    data = data[['iGC', 'jGC', 'MONTH', 'DAY', 'LON', 'LAT', 'OBS', 'MOD',
                 'PREC', 'ALBEDO_SWIR', 'BLENDED_ALBEDO']]

    # Calculate model - observation
    data['DIFF'] = data['MOD'] - data['OBS']

    # Print out some statistics
    cols = ['LAT', 'LON', 'MONTH', 'MOD']
    print('MODEL MAXIMUM : ', data['MOD'].max())
    print('MODEL MINIMUM : ', data['MOD'].min())
    print('TROPOMI MAXIMUM : ', data['OBS'].max())
    print('TROPOMI MINIMUM : ', data['OBS'].min())
    print('DIFFERENCE MAXIMUM : ', np.abs(data['DIFF']).max())
    print('DIFFERENCE MEAN : ', np.mean(data['DIFF']))
    print('DIFFERENCE STD : ', np.std(data['DIFF']))

    # Save the data out
    print(f'Saving data in {output_dir}/{settings.year}.pkl')
    gc.save_obj(data, join(output_dir, f'{settings.year}.pkl'))

else:
    ## ----------------------------------------- ##
    ## Load data for the year
    ## ----------------------------------------- ##
    print(f'Opening data in {output_dir}/{settings.year}.pkl')
    data = gc.load_obj(join(data_dir, prior_run))

print('Data is loaded.')
print(f'm = {data.shape[0]}')

## ----------------------------------------- ##
## Make and save observational mask
## ----------------------------------------- ##
# Create a vector for storing the observational filter
obs_filter = np.ones(data.shape[0], dtype=bool)
if filter_on_blended_albedo:
    BAF_filter = ((data['MONTH'].isin(np.arange(5, 11, 1))) |
                  (data['BLENDED_ALBEDO'] < blended_albedo_threshold))
    obs_filter = (obs_filter & BAF_filter)

if filter_on_albedo:
    albedo_filter = (data['ALBEDO_SWIR'] > albedo_threshold)
    obs_filter = (obs_filter & albedo_filter)

if filter_on_seasonal_latitude:
    latitude_filter = ((data['MONTH'].isin(np.arange(3, 12, 1))) |
                       (data['LAT'] <= 50))
    obs_filter = (obs_filter & latitude_filter)

obs_filter = pd.DataFrame({'MONTH' : data['MONTH'], 'FILTER' :  obs_filter})
obs_filter.to_csv(join(output_dir, 'obs_filter.csv'))

## ----------------------------------------- ##
## Additional data corrections
## ----------------------------------------- ##
if filter_on_blended_albedo:
    # Plot values
    if plot_dir is not None:
        fig, ax = fp.get_figax(rows=3, cols=4, aspect=1.25, sharex=True,
                               sharey=True)
        plt.subplots_adjust(wspace=0.1, hspace=0.7)
        for m in settings.months:
            d = data[data['MONTH'] == m]
            # fig, ax = fp.get_figax(aspect=1.75)
            axis = ax.flatten()[m-1]
            axis = fp.add_title(axis, cal.month_name[m])

            # Plot
            c = axis.hexbin(d['BLENDED_ALBEDO'], d['OBS'],
                          cmap=fp.cmap_trans('plasma_r'),
                          bins=np.arange(0, 2000),
                          vmin=0, vmax=2000)
            if m not in np.arange(5, 11, 1):
                axis.axvline(blended_albedo_threshold, c=fp.color(7), ls='--')

            # Axis labels
            if m in [1, 5]:
                axis = fp.add_labels(axis, '', 'XCH4\n(ppb)')
                axis.tick_params(labelbottom=False)
            elif m in [10, 11, 12]:
                axis = fp.add_labels(axis, 'Blended\nAlbedo', '')
                axis.tick_params(labelleft=False)
            elif m == 9:
                axis = fp.add_labels(axis, 'Blended\nAlbedo', 'XCH4\n(ppb)')
            else:
                axis.tick_params(labelbottom=False, labelleft=False)

            # Limits
            axis.set_xlim(0, 2)
            axis.set_ylim(1750, 1950)

        cax = fp.add_cax(fig, ax)
        cbar = fig.colorbar(c, cax=cax, ticks=500*np.arange(0, 5))
        cbar = fp.format_cbar(cbar, 'Count')
        # ax = fp.add_labels(ax, 'Blended Albedo', 'XCH4 (ppb)')

        fp.save_fig(fig, plot_dir,
                    f'blended_albedo_filter{suffix}')
        plt.close()


    # Apply the blended albedo filter
    BAF_filter = ((data['MONTH'].isin(np.arange(5, 11, 1))) |
                  (data['BLENDED_ALBEDO'] < blended_albedo_threshold))
    data = apply_filter(data, BAF_filter, 'blended albedo')

if filter_on_albedo:
    albedo_filter = (data['ALBEDO_SWIR'] > albedo_threshold)
    data = apply_filter(data, albedo_filter, 'albedo')

if filter_on_seasonal_latitude:
    latitude_filter = ((data['MONTH'].isin(np.arange(3, 12, 1))) |
                       (data['LAT'] <= 50))
    data = apply_filter(data, latitude_filter, 'winter latitude')

if remove_latitudinal_bias:
    # Correct the latitudinal bias
    coef = p.polyfit(data['LAT'], data['DIFF'], deg=1)
    bias_correction = p.polyval(data['LAT'], coef)
    data['MOD'] -= bias_correction
    data['DIFF'] -= bias_correction

    # Print information
    print(f'\nData has latitudinal bias removed.')
    print(f'    y = {coef[0]:.2f} + {coef[1]:.2f}x')
    print(f'    Model statistics:')
    summ = data.groupby('MONTH')[['MOD', 'DIFF']]
    summ = pd.concat([summ.min(), summ.max()], axis=1)
    print(summ)

# Save out result
if (filter_on_blended_albedo or filter_on_albedo or
    filter_on_seasonal_latitude or remove_latitudinal_bias):
    print(f'Saving data in {output_dir}/{settings.year}_corrected.pkl')
    gc.save_obj(data, join(output_dir, f'{settings.year}_corrected.pkl'))

nobs = data.shape[0]

## ------------------------------------------------------------------------ ##
## Analyze data
## ------------------------------------------------------------------------ ##
if analyze_biases:
    ## Spatial bias
    # Generate a grid on which to save out average difference
    lats, lons = gc.create_gc_grid(*settings.lats, settings.lat_delta,
                                   *settings.lons, settings.lon_delta,
                                   centers=False, return_xarray=False)

    # Save nearest latitude and longitude centers
    data['LAT_CENTER'] = lats[gc.nearest_loc(data['LAT'].values, lats)]
    data['LON_CENTER'] = lons[gc.nearest_loc(data['LON'].values, lons)]

    # Group on that grid
    s_b = data.groupby(['LAT_CENTER', 'LON_CENTER']).mean()['DIFF']
    s_b = s_b.to_xarray().rename({'LAT_CENTER' : 'lats',
                                 'LON_CENTER' : 'lons'})
    print('Spatial bias analyzed.')

    ## Latitudinal bias
    data['LAT_BIN'] = pd.cut(data['LAT'], lat_bins)
    l_b = gc.group_data(data, groupby=['LAT_BIN'])
    print('Latitudinal bias analyzed.')

    ## Seasonality
    m_b = gc.group_data(data, groupby=['MONTH'])
    print('Monthly bias analyzed.')

    ## Latitudinal bias and seasonality
    data.loc[:, 'SEASON'] = 'DJF'
    data.loc[data['MONTH'].isin([3, 4, 5]), 'SEASON'] = 'MAM'
    data.loc[data['MONTH'].isin([6, 7, 8]), 'SEASON'] = 'JJA'
    data.loc[data['MONTH'].isin([9, 10, 11]), 'SEASON'] = 'SON'
    lm_b = gc.group_data(data, groupby=['LAT_BIN', 'SEASON'])
    print('Seasonal latitudinal bias analyzed.')

    ## Albedo
    data['ALBEDO_BIN'] = pd.cut(data['ALBEDO_SWIR'], albedo_bins)
    a_b = gc.group_data(data, groupby=['ALBEDO_BIN'])
    print('Albedo bias analyzed.')

    ## -------------------------------------------------------------------- ##
    ## Plot data
    ## -------------------------------------------------------------------- ##
    if plot_dir is not None:
        ## ----------------------------------------- ##
        ## Scatter plot
        ## ----------------------------------------- ##
        fig, ax, c = gc.plot_comparison(data['OBS'].values, data['MOD'].values,
                                        lims=[1750, 1950], vmin=0, vmax=3e4,
                                        xlabel='Observation', ylabel='Model')
        ax.set_xticks(np.arange(1750, 2000, 100))
        ax.set_yticks(np.arange(1750, 2000, 100))
        fp.save_fig(fig, plot_dir, f'prior_bias{suffix}')
        plt.close()

        ## ----------------------------------------- ##
        ## Spatial bias
        ## ----------------------------------------- ##
        fig, ax = fp.get_figax(maps=True, lats=s_b.lats, lons=s_b.lons)
        c = s_b.plot(ax=ax, cmap='RdBu_r', vmin=-30, vmax=30,
                     add_colorbar=False)
        cax = fp.add_cax(fig, ax)
        cb = fig.colorbar(c, ax=ax, cax=cax)
        cb = fp.format_cbar(cb, 'Model - Observation')
        ax = fp.format_map(ax, s_b.lats, s_b.lons)
        ax = fp.add_title(ax, 'Spatial Bias in Prior Run')
        fp.save_fig(fig, plot_dir, f'prior_spatial_bias{suffix}')
        plt.close()

        ## ----------------------------------------- ##
        ## Latitudinal bias
        ## ----------------------------------------- ##
        l_b['LAT'] = l_b['LAT_BIN'].apply(lambda x: x.mid)
        fig, ax = fp.get_figax(aspect=1.75)
        ax.errorbar(l_b['LAT'], l_b['mean'], yerr=l_b['std'],
                    color=fp.color(4))
        ax.set_xticks(np.arange(10, 70, 10))
        ax.set_xlim(10, 60)
        ax = fp.add_labels(ax, 'Latitude', 'Model - Observation')
        ax = fp.add_title(ax, 'Latitudinal Bias in Prior Run')
        fp.save_fig(fig, plot_dir, f'prior_latitudinal_bias{suffix}')

        ## ----------------------------------------- ##
        ## Monthly bias
        ## ----------------------------------------- ##
        m_b['month'] = pd.to_datetime(m_b['MONTH'], format='%m')
        m_b['month'] = m_b['month'].dt.month_name().str[:3]
        fig, ax = fp.get_figax(aspect=1.75)
        ax.errorbar(m_b['month'], m_b['mean'], yerr=m_b['std'],
                    color=fp.color(4))
        ax = fp.add_labels(ax, 'Month', 'Model - Observation')
        ax = fp.add_title(ax, 'Seasonal Bias in Prior Run')
        fp.save_fig(fig, plot_dir, f'prior_seasonal_bias{suffix}')
        plt.close()

        ## ----------------------------------------- ##
        ## Seasonal latitudinal bias
        ## ----------------------------------------- ##
        lm_b['LAT'] = lm_b['LAT_BIN'].apply(lambda x: x.mid)
        fig, ax = fp.get_figax(aspect=1.75)
        ax.errorbar(l_b['LAT'], l_b['mean'], yerr=l_b['std'],
                    color=fp.color(4))
        ax.set_xticks(np.arange(10, 70, 10))
        ax.set_xlim(10, 60)
        ax = fp.add_labels(ax, 'Latitude', 'Model - Observation')
        ax = fp.add_title(ax, f'Latitudinal Bias in Prior Run')
        linestyles = ['solid', 'dotted', 'dashed', 'dashdot']
        for i, season in enumerate(np.unique(lm_b['SEASON'])):
            d = lm_b[lm_b['SEASON'] == season]
            ax.plot(d['LAT'].values, d['mean'].values, color=fp.color(4),
                    label=season, ls=linestyles[i], lw=0.5)
        fp.add_legend(ax)
        fp.save_fig(fig, plot_dir,
                    f'prior_seasonal_latitudinal_bias{suffix}')
        plt.close()

        ## ----------------------------------------- ##
        ## Albedo bias
        ## ----------------------------------------- ##
        a_b['ALBEDO'] = a_b['ALBEDO_BIN'].apply(lambda x: x.mid)
        fig, ax = fp.get_figax(aspect=1.75)
        ax.errorbar(a_b['ALBEDO'], a_b['mean'], yerr=a_b['std'],
                    color=fp.color(4))
        ax.set_xticks(np.arange(0, 1, 0.2))
        ax = fp.add_labels(ax, 'Albedo', 'Model - Observation')
        ax = fp.add_title(ax, 'Albedo Bias in Prior Run')
        fp.save_fig(fig, plot_dir, f'prior_albedo_bias{suffix}')
        plt.close()

## ------------------------------------------------------------------------ ##
## Calculate So
## ------------------------------------------------------------------------ ##
if calculate_so:
    # We calculate the mean bias, observation, and precision on the GEOS-Chem
    # grid, accounting for the squaring of the precision
    groupby = ['iGC', 'jGC', 'MONTH']
    group_quantities = ['DIFF', 'OBS', 'PREC_SQ']
    data['PREC_SQ'] = data['PREC']**2
    res_err = data.groupby(groupby).mean()[group_quantities].reset_index()
    res_err['PREC_SQ'] **= 0.5

    # Rename the columns
    res_err = res_err.rename(columns={'DIFF' : 'AVG_DIFF',
                                      'OBS' : 'AVG_OBS',
                                      'PREC_SQ' : 'AVG_PREC'})

    # Merge this data back into the original data frame
    data = pd.merge(data, res_err, on=groupby, how='left')

    # Subtract the bias from the difference to calculate the residual error
    # This is equivalent to ZQ's eps quantity
    data['RES_ERR'] = data['DIFF'] - data['AVG_DIFF']

    # Next we calculate the average residual error
    avg_err = data.groupby(groupby).mean()['RES_ERR'].reset_index()
    avg_err = avg_err.rename(columns={'RES_ERR' : 'AVG_RES_ERR'})
    # ZQ saves out:
    # average residual error as err_month.pkl,
    # average observations as obs_month.pkl
    # average precision as prec_month.pkl

    # Now calculate the gridded variance and standard deviation of the
    # residual error. The standard deviation is weighted by the number
    # of observations in a grid cell because this will decrease the
    # error in a grid cell.
    # (sigma_squared and rrsd, respectively, in ZQ's code)
    data = pd.merge(data, avg_err, on=groupby, how='left')
    data['VAR'] = (data['RES_ERR'] - data['AVG_RES_ERR'])**2
    var = data.groupby(groupby).mean()[['VAR', 'OBS']].reset_index()
    var = var.rename(columns={'OBS' : 'AVG_OBS'})
    var['STD'] = var['VAR']**0.5/var['AVG_OBS'] # rrsd
    # var['VAR'] is sigmasq

    # Merge these final variances back into the data (first removing
    # the initial variance calculation, since this was an intermediary)
    data = data.drop(columns=['VAR', 'AVG_OBS'])
    data = pd.merge(data, var, on=groupby, how='left')

    # Scale by the observations
    data['STD'] *= data['OBS']
    data['VAR'] = data['STD']**2

    # Where the variance calculated by the residual error method is less
    # than the precision squared value calculated above, set the error equal
    # to precision squared
    cond = data['VAR'] < data['PREC_SQ']
    print(f'We replace {cond.sum()} instances where the residual error is less than the instrumental error.')
    data.loc[:, 'SO'] = data['VAR']
    data.loc[cond, 'SO'] = data.loc[cond, 'PREC_SQ']

    # and then update std
    data.loc[:, 'STD'] = data['SO']**0.5

    err_mean = data['STD'].mean()
    print(f'We find a mean error of {err_mean:.2f} ppb.' )

    # Save out the data
    print(f'Saving data in {output_dir}/{settings.year}_corrected.pkl')
    gc.save_obj(data, join(output_dir, f'{settings.year}_corrected.pkl'))

## ------------------------------------------------------------------------ ##
## Plots
## ------------------------------------------------------------------------ ##
if (plot_dir is not None) and calculate_so:
    plot_data = copy.deepcopy(data[['iGC', 'jGC', 'MONTH', 'LON', 'LAT',
                                    'OBS', 'MOD', 'DIFF', 'PREC',
                                    'ALBEDO_SWIR', 'BLENDED_ALBEDO', 'SO']])

    # Generate a grid on which to save out average difference
    lats, lons = gc.create_gc_grid(*settings.lats, settings.lat_delta,
                                   *settings.lons, settings.lon_delta,
                                   centers=False, return_xarray=False)

    # Save nearest latitude and longitude centers
    lat_idx = gc.nearest_loc(plot_data['LAT'].values, lats)
    plot_data.loc[:, 'LAT_CENTER'] = lats[lat_idx]
    lon_idx = gc.nearest_loc(plot_data['LON'].values, lons)
    plot_data.loc[:, 'LON_CENTER'] = lons[lon_idx]

    # (and seasonally, just to show the variability in coverage)
    plot_data.loc[:, 'SEASON'] = 'DJF'
    plot_data.loc[plot_data['MONTH'].isin([3, 4, 5]), 'SEASON'] = 'MAM'
    plot_data.loc[plot_data['MONTH'].isin([6, 7, 8]), 'SEASON'] = 'JJA'
    plot_data.loc[plot_data['MONTH'].isin([9, 10, 11]), 'SEASON'] = 'SON'

    # Also calculate seasonal errors
    # We calculate the mean bias, observation, and precision on the GEOS-Chem
    # grid, accounting for the squaring of the precision
    groupby = ['LAT_CENTER', 'LON_CENTER', 'SEASON']
    group_quantities = ['DIFF', 'OBS', 'PREC_SQ']
    plot_data['PREC_SQ'] = plot_data['PREC']**2
    res_err = plot_data.groupby(groupby).mean()[group_quantities].reset_index()
    res_err['PREC_SQ'] **= 0.5

    # Rename the columns
    res_err = res_err.rename(columns={'DIFF' : 'AVG_DIFF',
                                      'OBS' : 'AVG_OBS',
                                      'PREC_SQ' : 'AVG_PREC'})

    # Merge this plot_data back into the original plot_data frame
    plot_data = pd.merge(plot_data, res_err, on=groupby, how='left')

    # Subtract the bias from the difference to calculate the residual error
    # This is equivalent to ZQ's eps quantity
    plot_data['RES_ERR'] = plot_data['DIFF'] - plot_data['AVG_DIFF']

    # Next we calculate the average residual error
    avg_err = plot_data.groupby(groupby).mean()['RES_ERR'].reset_index()
    avg_err = avg_err.rename(columns={'RES_ERR' : 'AVG_RES_ERR'})

    # Now calculate the gridded variance and standard deviation of the
    # residual error. The standard deviation is weighted by the number
    # of observations in a grid cell because this will decrease the
    # error in a grid cell.
    # (sigma_squared and rrsd, respectively, in ZQ's code)
    plot_data = pd.merge(plot_data, avg_err, on=groupby, how='left')
    plot_data['VAR'] = (plot_data['RES_ERR'] - plot_data['AVG_RES_ERR'])**2
    d_p = plot_data.groupby(groupby).agg({'VAR' : 'mean',
                                          'OBS' : 'mean',
                                          'DIFF' : 'count'})

    d_p = d_p.rename(columns={'OBS' : 'AVG_OBS', 'DIFF' : 'COUNT'})
    d_p['STD'] = d_p['VAR']**0.5#/d_p['AVG_OBS']
    d_p = d_p[['STD', 'AVG_OBS', 'COUNT']].to_xarray().rename({'LAT_CENTER' : 'lats',
                                                      'LON_CENTER' : 'lons'})

    fig, ax = fp.get_figax(rows=1, cols=4, maps=True,
                           lats=d_p.lats, lons=d_p.lons)
    fig_c, ax_c = fp.get_figax(rows=1, cols=4, maps=True,
                               lats=d_p.lats, lons=d_p.lons)
    fig_e, ax_e = fp.get_figax(rows=1, cols=4, maps=True,
                               lats=d_p.lats, lons=d_p.lons)
    for i, s in enumerate(['DJF', 'MAM', 'JJA', 'SON']):
        d = d_p.where(d_p.SEASON == s, drop=True)

        c = d['AVG_OBS'].plot(ax=ax[i], cmap='plasma', vmin=1800, vmax=1900,
                              add_colorbar=False)
        c_e = d['STD'].plot(ax=ax_e[i], cmap='plasma', vmin=0, vmax=25,
                            add_colorbar=False)
        c_c = d['COUNT'].plot(ax=ax_c[i], cmap='afmhot', vmin=0, vmax=300,
                              add_colorbar=False)
        for j, axis in enumerate([ax[i], ax_e[i], ax_c[i]]):
            axis = fp.format_map(axis, d.lats, d.lons)
            axis = fp.add_title(axis, s)
    cax = fp.add_cax(fig, ax)
    cb = fig.colorbar(c, ax=ax, cax=cax)
    cb = fp.format_cbar(cb, 'XCH4\n(ppb)')
    fp.save_fig(fig, plot_dir, f'observations{suffix}')

    cax_e = fp.add_cax(fig_e, ax_e)
    cb_e = fig.colorbar(c_e, ax=ax_e, cax=cax_e)
    cb_e = fp.format_cbar(cb_e, 'St. Dev.\n(ppb)')
    fp.save_fig(fig_e, plot_dir, f'errors{suffix}')

    cax_c = fp.add_cax(fig_c, ax_c)
    cb_c = fig.colorbar(c_c, ax=ax_c, cax=cax_c)
    cb_c = fp.format_cbar(cb_c, 'Count')
    fp.save_fig(fig_c, plot_dir, f'counts{suffix}')

    # Now plot the histogramS

    fig, ax = fp.get_figax(aspect=1.75)
    ax.hist(data['STD'], bins=250, density=True, color=fp.color(4))
    ax.set_xlim(0, 25)
    ax = fp.add_labels(ax, 'Observational Error (ppb)', 'Count')
    ax = fp.add_title(ax, 'Observational Error')
    fp.save_fig(fig, plot_dir, 'observational_error')

    # SEASONAL
    fig, ax = fp.get_figax(aspect=1.75)
    ax.set_xlim(0, 25)
    ax = fp.add_labels(ax, 'Observational Error (ppb)', 'Count')
    ax = fp.add_title(ax, 'Observational Error')
    for i, season in enumerate(np.unique(data['SEASON'])):
        hist_data = data[data['SEASON'] == season]['STD']
        ax.hist(hist_data, histtype='step', bins=275, label=season,
                color=fp.color(2+2*i), lw=1)
        ax.axvline(hist_data.mean(), color=fp.color(2+2*i), lw=1, ls=':')
    ax = fp.add_legend(ax)
    fp.save_fig(fig, plot_dir, 'observational_error_seasonal')

    # LATITUDE
    fig, ax = fp.get_figax(aspect=1.75)
    ax.set_xlim(0, 25)
    ax = fp.add_labels(ax, 'Observational Error (ppb)', 'Count')
    ax = fp.add_title(ax, 'Observational Error')
    for i, lat_bin in enumerate(np.unique(data['LAT_BIN'])):
        hist_data = data[data['LAT_BIN'] == lat_bin]['STD']
        ax.hist(hist_data, histtype='step', bins=275, label=lat_bin,
                color=fp.color(2*i), lw=1)
        ax.axvline(hist_data.mean(), color=fp.color(2*i), lw=1, ls=':')

    ax = fp.add_legend(ax)
    fp.save_fig(fig, plot_dir, 'observational_error_latitude_hist')

    fig, ax = fp.get_figax(aspect=1.75)
    ax.scatter(data['LAT'], data['STD'], c=fp.color(4), s=2, alpha=0.1)
    ax = fp.add_labels(ax, 'Latitude', 'Observational Error (ppb)')
    ax = fp.add_title(ax, 'Observational Error')
    fp.save_fig(fig, plot_dir, 'observational_error_latitude_scatter')

    # ALBEDO
    fig, ax = fp.get_figax(aspect=1.75)
    ax.set_xlim(0, 25)
    ax = fp.add_labels(ax, 'Observational Error (ppb)', 'Count')
    ax = fp.add_title(ax, 'Observational Error')
    for i, alb_bin in enumerate(np.unique(data['ALBEDO_BIN'])):
        hist_data = data[data['ALBEDO_BIN'] == alb_bin]['STD']
        ax.hist(hist_data, histtype='step', bins=275, label=alb_bin,
                color=fp.color(2*i), lw=1)
        ax.axvline(hist_data.mean(), color=fp.color(2*i), lw=1, ls=':')

    ax = fp.add_legend(ax)
    fp.save_fig(fig, plot_dir, 'observational_error_albedo_hist')

    fig, ax = fp.get_figax(aspect=1.75)
    ax.scatter(data['ALBEDO_SWIR'], data['STD'], c=fp.color(4), s=2, alpha=0.1)
    ax = fp.add_labels(ax, 'Albedo', 'Observational Error (ppb)')
    ax = fp.add_title(ax, 'Observational Error')
    fp.save_fig(fig, plot_dir, 'observational_error_albedo_scatter')

## ------------------------------------------------------------------------ ##
## Save out inversion quantities
## ------------------------------------------------------------------------ ##
print(f'Saving data in {output_dir}')

y = xr.DataArray(data['OBS'], dims=('nobs'))
y.to_netcdf(join(output_dir, 'y.nc'))

ya = xr.DataArray(data['MOD'], dims=('nobs'))
ya.to_netcdf(join(output_dir, 'ya.nc'))

if calculate_so:
    so = xr.DataArray(data['SO'], dims=('nobs'))
    so.to_netcdf(join(output_dir, 'so.nc'))

print(f'Number of observations : {nobs}')
print('=== CODE COMPLETE ====')
