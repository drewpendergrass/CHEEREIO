import sys
import pytest
import xarray as xr
import numpy as np
import json
sys.path.append('../core/')
import observation_operators as obsop
import tropomi_tools as tt
import testing_tools

#Test that the methane qa values are all high quality
def test_all_good_qa_values_in_methane():
	met = tt.read_tropomi('data_for_tests/METHANE_TEST/tropomi/S5P_OFFL_L2__CH4____20190107T062928_20190107T081058_06397_01_010202_20190113T084352.nc', 'CH4', filterinfo=None, calcError = False)
	assert np.sum(met['qa_value']>0.5) == len(met['qa_value'])
