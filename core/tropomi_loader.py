import tropomi_tools as tt
from datetime import datetime

trans = tt.TROPOMI_Translator()
trans.initialReadDate()
timeperiod = [datetime.strptime('2019-01-01 00:00', '%Y-%m-%d %H:%M'), datetime.strptime('2019-01-07 23:59', '%Y-%m-%d %H:%M')]

trop_dat = trans.getTROPOMI('CH4',timeperiod)