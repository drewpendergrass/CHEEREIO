import sys
import pytest
sys.path.append('../')
from GC_Translator import GC_Translator

#These tests ensure that we are subsetting columns correctly in the GC_Translator class.

def test_col_subset():
	gt = GC_Translator('data_for_tests/TEST_0001/','20190101_0000',computeStateVec = True)
	