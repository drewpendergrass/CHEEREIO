from glob import glob
import toolbox as tx 
import sys 

cmdarg = str(sys.argv[1])
if cmdarg=="TESTING":
	data = tx.getSpeciesConfig(testing=True)
	latgrid,longrid = tx.getLatLonVals(data,True)
else:
	data = tx.getSpeciesConfig(testing=False)
	latgrid,longrid = tx.getLatLonVals(data,False)

path_to_scratch = f"{data['MY_PATH']}/{data['RUN_NAME']}/scratch"
columns = glob(f'{path_to_scratch}/**/*.npy',recursive=True)
numcols = len(columns)
num_cells = len(latgrid)*len(longrid)

if numcols==num_cells:
	with open(f"{path_to_scratch}/ALL_COLUMNS_FOUND", "w") as f:
		f.write("All columns found.\n") #If so, save flag file to ensemble folder
