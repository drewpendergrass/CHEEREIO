import pytest
import sys
sys.path.append('../core/')
import testing_tools

def pytest_sessionfinish(session, exitstatus):
	#turn off override so that CHEEREIO continues to behave as expected.
	testing_tools.turnOffOverride()