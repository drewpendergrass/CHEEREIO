import sys
import pytest
import xarray as xr
import numpy as np
sys.path.append('../')
from GC_Translator import GC_Translator

#These tests ensure that we are subsetting columns correctly in the GC_Translator class.

def test_col_subset():
	gt = GC_Translator('data_for_tests/TEST_0001/','20190101_0000',computeStateVec = True)
	ds = xr.load_dataset('data_for_tests/TEST_0001/GEOSChem.Restart.20190101_0000z.nc4')
	da = np.array(ds[f'SpeciesRst_CH4']).squeeze()
