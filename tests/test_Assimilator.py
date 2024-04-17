import sys
import pytest
import xarray as xr
import numpy as np
import json
import scipy.linalg as la
sys.path.append('../core/')
from Assimilator import Assimilator
import testing_tools

#Just test the LETKF math -- don't bother with making sure everything else loads right
#This test also implicitly makes sure that assimilator constructor works.
def test_LETKF_calculation():
	testing_tools.setupPytestSettings('methane')
	assim = testing_tools.prepTestAssimilator()
	assim_answer = assim.analysisEnsemble
	#calculate actual answer, using assumption of 2 ensemble members and 0 inflation.
	ptilde = la.inv(np.diag(np.ones(2)) + (np.transpose(assim.Ypert_background)@la.inv(assim.R)@assim.Ypert_background))
	bigW = la.sqrtm(ptilde)
	lilW = ptilde @ np.transpose(assim.Ypert_background)@la.inv(assim.R)@assim.ydiff
	real_answer = np.zeros(np.shape(assim.Xpert_background))
	for i in range(np.shape(assim_answer)[1]):
		w = bigW[:,i]+lilW
		real_answer[:,i] = (assim.Xpert_background@w) + assim.xbar_background
	assert np.allclose(assim_answer,real_answer)

#Check that localized statevector puts emission sf at end of column
def testStateVecSF():
	testing_tools.setupPytestSettings('methane')
	assim = Assimilator('20190108_0000',1,1)
	sv = assim.combineEnsemble(19,43) #StateVec for a random location
	colinds = assim.gt[1].getColumnIndicesFromLocalizedStateVector(19,43) #Column indices of localized statevec
	sf_from_statevec = sv[colinds,0][-1] #last column entry in localized statevec should be the emission sf for ensemble member 1
	assert np.abs(sf_from_statevec-assim.gt[1].getEmisSF('CH4')[19,43])<1e-16 #check the above

#Test that we are correctly scaling emissions to match initial standard deviation, if this is the only active analysis correction
def test_LETKF_emis_SF_scaling():
	testing_tools.setupPytestSettings('methane')
	assim = testing_tools.prepTestAssimilator(59,101)
	errors = []
	#Test that we successfully update SF if std collapses
	analysisSubset,backgroundSubset = assim.getAnalysisAndBackgroundColumn(59,101,doBackground=True,doPerts=False) #Get column subsets
	analysisSubset[-1,:] = np.array([1.2,1.3]) #Set posterior scale factors
	backgroundSubset[-1,:] = np.array([1.0,1.5]) #Set prior scale factors
	assim.InitEmisSTD['CH4'][59,101] = 0.6 #assume initial background std was 0.6
	assim = testing_tools.setupAssimilatorForAnalysisCorrectionUnitTest(assim,'InflateScalingsToXOfInitialStandardDeviation',0.3)
	correctedAnalysisSubset = assim.applyAnalysisCorrections(analysisSubset,backgroundSubset,59,101)
	trueAnalysisSubset = analysisSubset
	trueAnalysisSubset[-1,:] = ((np.array([1.2,1.3])-1.25)*3.6) + 1.25 
	if not np.allclose(trueAnalysisSubset,correctedAnalysisSubset):
		errors.append('Assimilator failed to update analysis scale factors to X of initial standard deviation')
	#Test that no update happens if std does not collapse
	analysisSubset,backgroundSubset = assim.getAnalysisAndBackgroundColumn(59,101,doBackground=True,doPerts=False) #Get column subsets
	analysisSubset[-1,:] = np.array([1.2,1.3]) #Set posterior scale factors
	backgroundSubset[-1,:] = np.array([1.0,1.5]) #Set prior scale factors
	assim.InitEmisSTD['CH4'][59,101] = 0.1 #assume initial background std was 0.1
	assim = testing_tools.setupAssimilatorForAnalysisCorrectionUnitTest(assim,'InflateScalingsToXOfInitialStandardDeviation',0.3)
	correctedAnalysisSubset = assim.applyAnalysisCorrections(analysisSubset,backgroundSubset,59,101)
	if not np.allclose(analysisSubset,correctedAnalysisSubset):
		errors.append('Assimilator failed to leave scale factors at analysis standard deviation when this behavior was expected.')
	assert not errors, "errors occured:\n{}".format("\n".join(errors))     


#Test that we are correctly scaling emissions to match initial standard deviation, if this is the only active analysis correction
def test_MinMaxSF():
	testing_tools.setupPytestSettings('methane')
	assim = testing_tools.prepTestAssimilator(59,101)
	errors = []
	#Test that we successfully update SF if std collapses
	analysisSubset,backgroundSubset = assim.getAnalysisAndBackgroundColumn(59,101,doBackground=True,doPerts=False) #Get column subsets
	analysisSubset[-1,:] = np.array([-0.3,0.2]) #Set posterior scale factors
	backgroundSubset[-1,:] = np.array([0.1,0.3]) #Set prior scale factors
	assim = testing_tools.setupAssimilatorForAnalysisCorrectionUnitTest(assim,'MinimumScalingFactorAllowed',0.01)
	correctedAnalysisSubset = assim.applyAnalysisCorrections(analysisSubset,backgroundSubset,59,101)
	trueAnalysisSubset = analysisSubset
	trueAnalysisSubset[-1,:] = np.array([0.01,0.2])
	if not np.allclose(trueAnalysisSubset,correctedAnalysisSubset):
		errors.append('Assimilator failed to enforce minimum scale factor correctly.')
	#Test that no update happens if std does not collapse
	analysisSubset,backgroundSubset = assim.getAnalysisAndBackgroundColumn(59,101,doBackground=True,doPerts=False) #Get column subsets
	analysisSubset[-1,:] = np.array([12.2,25.3]) #Set posterior scale factors
	backgroundSubset[-1,:] = np.array([19.2,15.6]) #Set prior scale factors
	assim = testing_tools.setupAssimilatorForAnalysisCorrectionUnitTest(assim,'MaximumScalingFactorAllowed',20)
	correctedAnalysisSubset = assim.applyAnalysisCorrections(analysisSubset,backgroundSubset,59,101)
	trueAnalysisSubset = analysisSubset
	trueAnalysisSubset[-1,:] = np.array([12.2,20])
	if not np.allclose(trueAnalysisSubset,correctedAnalysisSubset):
		errors.append('Assimilator failed to enforce maximum scale factor correctly.')
	assert not errors, "errors occured:\n{}".format("\n".join(errors))     


