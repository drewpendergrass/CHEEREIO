from glob import glob
import toolbox as tx 
import sys 

data = tx.getSpeciesConfig()
latgrid,longrid = tx.getLatLonVals(data)

path_to_scratch = f"{data['MY_PATH']}/{data['RUN_NAME']}/scratch"
columns = glob(f'{path_to_scratch}/**/*.npy',recursive=True)
numcols = len(columns)
num_cells = len(latgrid)*len(longrid)

if numcols==num_cells:
	with open(f"{path_to_scratch}/ALL_COLUMNS_FOUND", "w") as f:
		f.write("All columns found.\n") #If so, save flag file to ensemble folder
