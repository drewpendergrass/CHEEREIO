import sys
import pytest
import xarray as xr
import numpy as np
import json
sys.path.append('../core/')
from Assimilator import Assimilator
import testing_tools

def test_assim_loads():
	#override ens_config so that we are set to interpret the data properly
	testing_tools.setupPytestSettings('methane')
	assim = Assimilator('20190108_0000',1,1)
	assert 1 == 1

