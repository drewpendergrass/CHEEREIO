import tropomi_tools as tt
from datetime import datetime
from glob import glob
import pickle

trans = tt.TROPOMI_Translator()
timeperiod = [datetime.strptime('2019-01-01 00:00', '%Y-%m-%d %H:%M'), datetime.strptime('2019-01-07 23:59', '%Y-%m-%d %H:%M')]

trop_dat = trans.getTROPOMI('CH4',timeperiod)

# sourcedir = trans.spc_config['TROPOMI_dirs']['CH4']
# obs_list = glob(f'{sourcedir}/S5P_*.nc')
# obs_list.sort()

# with open(f"{trans.scratch}/tropomi_dates.pickle", 'rb') as handle:
# 	TROPOMI_date_dict = pickle.load(handle)

# obs_dates = TROPOMI_date_dict['CH4']
