import numpy as np
from glob import glob
import sys

#Calculates the time complexity of the CHEEREIO run, by parsing the log files.
#USAGE: python grab_time_info_from_logs.py PATH/TO/LOG/DIRECTORY
log_directory = str(sys.argv[1])

def getFile(filename):
	with open(filename) as f:
		lines = f.readlines()
	return lines

def grabGCWallTime(lines):
	times = []
	for line in lines:
		if line.startswith('real'):
			times.append(float(line.split()[1]))
	times = np.array(times)
	return times

def grabLETKFWallTime(lines):
	ens_gather_times = []
	ens_save_times = []
	isLETKF = False #are we in LETKF mode, as opposed to scale mode
	for line in lines:
		if line.startswith('Core'):
			splitline = line.split()
			if splitline[2] == "gathered":
				ens_gather_times.append(float(splitline[5]))
				if splitline[8] == "LETKF":
					isLETKF = True
				elif splitline[8] == "scaling":
					isLETKF = False
			if (splitline[2] == "completed") and isLETKF:
				ens_save_times.append(float(splitline[10]))
	ens_gather_times = np.array(ens_gather_times)
	ens_save_times = np.array(ens_save_times)
	return [ens_gather_times,ens_save_times]


def grabLETKFMasterTime(lines):
	ens_gather_times = []
	for line in lines:
		if line.startswith('Core gathered'):
			ens_gather_times.append(float(line.split()[6]))
	ens_gather_times = np.array(ens_gather_times)
	return ens_gather_times

def meanGCWallTime(log_directory):
	gc_log_files = glob(f'{log_directory}/ensemble_slurm_*.err')
	gc_log_files.sort()
	times_list = []
	for file in gc_log_files:
		lines = getFile(file)
		times_list.append(grabGCWallTime(lines))
	times = np.concatenate(times_list)
	mean_time = np.mean(times)
	return mean_time

def meanLETKFTime(log_directory):
	letkf_files = glob(f'{log_directory}/letkf_*.out')
	letkf_files.sort()
	ens_gather_times_list = []
	ens_save_times_list = []
	for file in letkf_files:
		lines = getFile(file)
		ens_gather_times,ens_save_times = grabLETKFWallTime(lines)
		ens_gather_times_list.append(ens_gather_times)
		ens_save_times_list.append(ens_save_times)
	egt = np.concatenate(ens_gather_times_list)
	est = np.concatenate(ens_save_times_list)
	mean_gather_time = np.mean(egt)
	mean_save_time = np.mean(est)
	return [mean_gather_time,mean_save_time]

def meanLETKFMasterTime(log_directory):
	file = f'{log_directory}/letkf_master.out'
	lines = getFile(file)
	times = grabLETKFMasterTime(lines)
	mean_time = np.mean(times)
	return mean_time

gc_wall_time = meanGCWallTime(log_directory)
letkf_gather_time,letkf_compute_time = meanLETKFTime(log_directory)
letkf_master_gather_time = meanLETKFMasterTime(log_directory)
print(f'Mean GC wall time: {np.round(gc_wall_time/60,2)} minutes.')
print(f'Mean time to gather ensemble for LETKF: {np.round(letkf_gather_time/60,2)} minutes.')
print(f'Mean time to compute LETKF and save columns: {np.round(letkf_compute_time/60,2)} minutes.')
print(f'Mean time to gather LETKF and columns to overwrite for next run: {np.round(letkf_master_gather_time/60,2)} minutes.')
print('')
total_time = gc_wall_time+letkf_gather_time+letkf_compute_time+letkf_master_gather_time
gc_frac = gc_wall_time/total_time
cheereio_frac = (total_time-gc_wall_time)/total_time
print(f'Percent of time spent in GEOS-Chem: {np.round(gc_frac,3)*100}%')
print(f'Percent of time spent in CHEEREIO: {np.round(cheereio_frac,3)*100}%')