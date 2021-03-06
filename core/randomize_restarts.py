from GC_Translator import GC_Translator
from glob import glob
import sys

#Utility script for testing. Randomly perturbs restart files in ensemble run directories
#While adding a bias factor to the "nature run."
print('Begin restart adjustment.')
path_to_ensemble = str(sys.argv[1]) #Emsemble directory; "/n/holyscratch01/jacob_lab/dpendergrass/GC-LETKF/Korea_NO2_test/ensemble_runs/"
timestamp = str(sys.argv[2]) #Restart timestamp; "20190101_0000"
biasshift = float(sys.argv[3]) #How much to change the bias.

subdirs = glob(f"{path_to_ensemble}/*/")
subdirs.remove(f"{path_to_ensemble}/logs/")
dirnames = [d.split('/')[-2] for d in subdirs]
subdir_numbers = [int(n.split('_')[-1]) for n in dirnames]
maxdir = max(subdir_numbers)
meanval = (maxdir+1)/2
for ens, directory in zip(subdir_numbers,subdirs):
	if ens!=0:
		print(f'Perturbing directory #{ens}.')
		gt = GC_Translator(directory, timestamp, False, True)
		gt.randomizeRestart(perturbation=0,bias=(ens/meanval)-1+biasshift)
		gt.saveRestart()

print('Restart adjustment complete.')
