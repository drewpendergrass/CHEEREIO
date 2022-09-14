import sys
import pytest
import xarray as xr
import numpy as np
import pandas as pd
import json
sys.path.append('../core/')
import observation_operators as obsop
import tropomi_tools as tt
import testing_tools

#Test that the methane qa values are all high quality
def test_all_good_qa_values_in_methane():
	met = tt.read_tropomi('data_for_tests/METHANE_TEST/tropomi/S5P_OFFL_L2__CH4____20190107T062928_20190107T081058_06397_01_010202_20190113T084352.nc', 'CH4', filterinfo=None, includeObsError = False)
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

#Make sure we average observations to GEOS-Chem grid correctly! Check that all 6 main parameters are averaged as expected.
def testAverageByGC():
	GC = testing_tools.makeMiniFakeDataSet()
	OBSDATA = testing_tools.makeMiniFakeObsData(latlocs = [4],lonlocs = [9],ntime = 10)
	iGC, jGC, tGC = obsop.nearest_loc(GC,OBSDATA)
	obsvals = np.arange(10)
	GCmappedtoobs = np.array([0,0,0,1,1,1,1,2,2,2])
	obsdata_results = obsop.averageByGC(iGC, jGC, tGC, GC,GCmappedtoobs,obsvals,doSuperObs=False)
	correctlon = np.repeat(10,3) #Nearest longitude to 9 is 10 on our 0,5,10 gc grid
	correctlat = np.repeat(5,3) #Nearest latitude to 4 is 5.
	testlat,testlon = obsdata_results.getLatLon()
	lonresult = np.array_equal(correctlon,testlon)
	latresult = np.array_equal(correctlat,testlat)
	testgc, testobs = obsdata_results.getCols()
	correctgc = np.arange(3)
	correctobs = np.array([np.mean(obsvals[0:3]), np.mean(obsvals[3:7]), np.mean(obsvals[7:10])])
	gcresult = np.array_equal(correctgc,testgc)
	obsresult = np.array_equal(correctobs,testobs)
	testtime = obsdata_results.getTime()
	correcttime = pd.date_range(start='2022-08-01',end='2022-08-08',periods=3).values.astype(int)
	timeresult = np.array_equal(testtime,correcttime)
	testnum = obsdata_results.getDataByKey('num_av')
	correctnum = np.array([3,4,3]) #3 matches, 4 matches, 3 matches
	numresult = np.array_equal(testnum,correctnum)
	assert lonresult and latresult and gcresult and obsresult and timeresult and numresult

def testSuperObsFunctions():
	GC = testing_tools.makeMiniFakeDataSet()
	OBSDATA = testing_tools.makeMiniFakeObsData(latlocs = [4],lonlocs = [9],ntime = 10)
	iGC, jGC, tGC = obsop.nearest_loc(GC,OBSDATA)
	obsvals = np.arange(10)
	obsInstrumentError = np.ones(10)*5
	modelTransportError = 5
	errorCorr = 0.1
	GCmappedtoobs = np.array([0,0,0,1,1,1,1,2,2,2])
	testobserr_default = obsop.averageByGC(iGC, jGC, tGC, GC,GCmappedtoobs,obsvals,doSuperObs=True,superObsFunction='default',obsInstrumentError = obsInstrumentError, modelTransportError = modelTransportError, errorCorr = errorCorr,minError=0).getDataByKey('err_av')
	correctobserr_default = [(np.mean(obsInstrumentError[0:3]) * np.sqrt(((1-errorCorr)/3) + errorCorr) )+modelTransportError,(np.mean(obsInstrumentError[3:7]) * np.sqrt(((1-errorCorr)/4) + errorCorr) )+modelTransportError,(np.mean(obsInstrumentError[7:10]) * np.sqrt(((1-errorCorr)/3) + errorCorr) )+modelTransportError]
	print(testobserr_default)
	print(correctobserr_default)
	obsresult_default = np.array_equal(testobserr_default,correctobserr_default)
	assert obsresult_default



