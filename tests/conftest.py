import pytest
import sys
import os
sys.path.append('../core/')
import testing_tools

def pytest_sessionfinish(session, exitstatus):
	#turn off override so that CHEEREIO continues to behave as expected.
	testing_tools.turnOffOverride()
	#Delete tropomi_dates file if it was produced.
	tropomi_dates_path = 'data_for_tests/METHANE_TEST/scratch/tropomi_dates.pickle'
	if os.path.exists(tropomi_dates_path):
		os.remove(tropomi_dates_path)

