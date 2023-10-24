import numpy as np
import xarray as xr
from animation_tools import *
import argparse

def str2bool(v):
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser(description='Animate a NetCDF file produced by CHEEREIO.')
parser.add_argument('-i', '--file_in', type=str, help='NetCDF file input (full path).')
parser.add_argument('-o', '--file_out', type=str, help='Movie file output (full path).')
parser.add_argument('-var', '--variable', type=str, help='Variable from NetCDF to plot.')
parser.add_argument('-log', '--lognormal_errors', type=str2bool, default=True, help='Lognormal errors?')
parser.add_argument('-bwr_cmap', '--bwr_cmap', type=str2bool, default=False, help='Use the blue white red colormap, centered on 1?')
parser.add_argument('-func', '--statistical_function', type=str, default='mean', help='Statistics to plot? mean, sd, max, min, or range.')
parser.add_argument('-anim_fps', '--anim_fps', type=int, default='8', help='Frames per second of animation.')
parser.add_argument('-cmin', '--cmin', type=float, default=np.nan, help='Minimum of colorbar? If nan (default) use data to infer.')
parser.add_argument('-cmax', '--cmax', type=float, default=np.nan, help='Maximum of colorbar? If nan (default) use data to infer.')
parser.add_argument('-start', '--start_date', type=str, default='None', help='Start date for animation? YYYYMMDD form, or "None" to ignore. Both start and end must be supplied if using.')
parser.add_argument('-end', '--end_date', type=str, default='None', help='End date for animation? YYYYMMDD form, or "None" to ignore')

# create dictionary of arguments
args = parser.parse_args()

if np.isnan(args.cmin):
	cmin = None
else:
	cmin = args.cmin

if np.isnan(args.cmax):
	cmax = None
else:
	cmax = args.cmax


ds = xr.open_dataset(args.file_in)

if args.start_date != 'None':
	start_date=datetime.strptime(args.start_date, "%Y%m%d")
	end_date=datetime.strptime(args.end_date, "%Y%m%d")
	ds = ds.sel(time=slice(start_date, end_date))

time = np.array(ds['time'])
timestr = [str(t)[0:16] for t in time]
lat = np.array(ds['lat'])
lon = np.array(ds['lon'])
ensmean = getEnsMean(func,variable,da,args.lognormal_errors)


animateData(m,ensmean,args.file_out,lon,lat,anim_fps = args.anim_fps, variable = args.variable,timestr = timestr, cmin = cmin, cmax=cmax, bwr_cmap=args.bwr_cmap)


