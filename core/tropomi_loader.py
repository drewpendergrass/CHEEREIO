import tropomi_tools as tt
from HIST_Translator import HIST_Translator
import xarray as xr

met = tt.read_tropomi('testNO2','NO2')
ht = HIST_Translator
GC = 