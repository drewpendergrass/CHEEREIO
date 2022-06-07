import sys
import pytest
import xarray as xr
import numpy as np
sys.path.append('../')
from GC_Translator import GC_Translator

#These tests ensure that we are subsetting columns correctly in the GC_Translator class.

#From the methane restart, localize the state vector about 10,10 then get the column from the localization.
#Compare this to the true column from the file
def test_col_subset_of_localized_state_vector_methane():
	gt = GC_Translator('data_for_tests/METHANE_TEST/TEST_0001/','20190101_0000',computeStateVec = True)
	locstatevec = gt.getStateVector(10,10)
	columninds = gt.getColumnIndicesFromLocalizedStateVector(10,10)
	column_from_statevec = locstatevec[columninds]
	ds = xr.load_dataset('data_for_tests/METHANE_TEST/TEST_0001/GEOSChem.Restart.20190101_0000z.nc4')
	da = np.array(ds[f'SpeciesRst_CH4']).squeeze()
	column_from_file = da[:,10,10]
	assert np.allclose(column_from_statevec,column_from_file)
