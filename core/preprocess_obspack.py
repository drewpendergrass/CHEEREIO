from obspack_tools import prep_obspack
import settings_interface as si 
import argparse
from datetime import datetime

spc_config = si.getSpeciesConfig()


parser = argparse.ArgumentParser(description='Preprocess obspack for GEOS-Chem and CHEEREIO; defaults to ens_config settings, but you can override.')
parser.add_argument('-i', '--raw_obspack_dir', type=str, default=spc_config['raw_obspack_path'], help='Input directory containing raw obspack data? Default: raw_obspack_path in ens_config.json')
parser.add_argument('-o', '--gc_obspack_dir', type=str, default=spc_config['gc_obspack_path'], help='Output directory for GC-formatted obspack data? Default: gc_obspack_path in ens_config.json')
parser.add_argument('-name', '--filename_format', type=str, default=spc_config['obspack_gc_input_file'], help='Output filename format for GC-formatted obspack data? Default: obspack_gc_input_file in ens_config.json')
parser.add_argument('-start', '--start_date', type=str, default=spc_config['START_DATE'], help='Start date for data processing (YYYYMMDD)? Default: START_DATE in ens_config.json')
parser.add_argument('-end', '--end_date', type=str, default=spc_config['END_DATE'], help='End date for data processing (YYYYMMDD)? Default: END_DATE in ens_config.json')

# create dictionary of arguments
args = parser.parse_args()
# Parse date
start_date = datetime.strptime(args.start_date, "%Y%m%d")
end_date = datetime.strptime(args.end_date, "%Y%m%d")

#Prep obspack
print(f'Loading obspack data from {args.raw_obspack_dir}...')
prep_obspack(args.raw_obspack_dir,args.gc_obspack_dir,args.filename_format,start_date,end_date)
print(f'Done! Parsed obspack data saved at {args.gc_obspack_dir}.')

