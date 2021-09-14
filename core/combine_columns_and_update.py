from glob import glob
import toolbox as tx 
import sys 

cmdarg = str(sys.argv[1])
if cmdarg=="TESTING":
	data = tx.getSpeciesConfig(testing=True)
else:
	data = tx.getSpeciesConfig(testing=False)
