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

#Test that we are correctly
def test_LETKF_emis_SF():
	testing_tools.setupPytestSettings('methane')
	assim = testing_tools.prepTestAssimilator(59,101)
	analysisSubset,backgroundSubset = assim.getAnalysisAndBackgroundColumn(59,101,doBackground=True,doPerts=False)
	assert True
	# colinds = assim.gt[1].getColumnIndicesFromLocalizedStateVector(59,101)

