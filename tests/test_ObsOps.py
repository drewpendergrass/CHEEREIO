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

#Check whether the indexing GC by observations works
def testGCindexing():
	GC = testing_tools.makeMiniFakeDataSet()
	OBSDATA = testing_tools.makeMiniFakeObsData(latlocs = [4],lonlocs = [9],ntime = 10)
	iGC, jGC, tGC = obsop.nearest_loc(GC,OBSDATA)
	#i is lon, j is lat, t is time. Expect all lats to be nearest middle index (1), all lons to be near the last index (2),
	#and for the first three and last three times to be on the edges (0 and 2) and the middle four to be in the middle (1)
	assert np.array_equal(jGC,np.repeat(1,10)) and np.array_equal(iGC,np.repeat(2,10)) and np.array_equal(tGC,np.array([0,0,0,1,1,1,1,2,2,2]))

#test that we get the correct GC columns for a given set of observations
def testGetGCCols():
	GC = testing_tools.makeMiniFakeDataSet()
	OBSDATA = testing_tools.makeMiniFakeObsData(latlocs = [4],lonlocs = [9],ntime = 10)
	results,_,_ = obsop.getGCCols(GC,OBSDATA,'TEST')
	results_level0 = results[:,0]
	#on the lat/lon grid at level 0, time zero, we get 5 (latind = 1 ==> second row. lonind = 2 ==> third column). 
	#Advance by 27 per timestep.
	correct_answer = 5+(27*np.array([0,0,0,1,1,1,1,2,2,2]))
	assert np.array_equal(results_level0,correct_answer)

