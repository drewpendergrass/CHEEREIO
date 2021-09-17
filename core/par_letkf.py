import letkf_utils as lu
import sys
import time

timestamp = str(sys.argv[1]) #Time to assimilate. Expected in form YYYYMMDD_HHMM, UTC time.
ensnum = int(sys.argv[2])
corenum = int(sys.argv[3])
testing = bool(sys.argv[4])

dateval = timestamp[0:4]+'-'+timestamp[4:6]+'-'+timestamp[6:8]

print(f'Core ({ensnum},{corenum}) is gathering ensemble at time {dateval}.')
start = time.time()
assimilator = lu.Assimilator(timestamp,ensnum,corenum,testing)
end = time.time()
print(f'Core ({ensnum},{corenum}) gathered ensemble in {end - start} seconds. Begin LETKF procedure.')
start = time.time()
assimilator.LETKF()
end = time.time()
print(f'Core ({ensnum},{corenum}) completed computation for {dateval} and saved columns in {end - start} seconds.')
# start = time.time()
# assimilator.updateRestartsAndScalingFactors()
# assimilator.saveRestartsAndScalingFactors()
# end = time.time()
# print(f'Saving completed in {end - start} seconds.')
# print(f'Core ({ensnum},{corenum}) processed its subset of data for {dateval}.')
print('-------------------END LETKF-------------------')

#THESE LINES ARE ONLY USED FOR TESTING
# import letkf_utils as lu
# import toolbox as tx
# import numpy as np
# timestamp = '20190101_0000'
# ensnum = 1
# corenum = 1
# testing = True
# assimilator = lu.Assimilator(timestamp,ensnum,corenum,testing)
# assimilator.prepareMeansAndPerts(23,0)
# assimilator.ybar_background-assimilator.Ypert_background[:,0]
# latinds,loninds = tx.getIndsOfInterest(23,0,testing=True)
# conc4D=assimilator.combineEnsembleForSpecies('NO2')
# obsEns = np.zeros([len(assimilator.ObsOp['NO2_SATELLITE'].obsinfo.getObsVal(23,0)),np.shape(conc4D)[3]])

# for i in range(16):
# 	obsEns[:,i],_,_ = assimilator.ObsOp['NO2_SATELLITE'].H(conc4D[:,:,:,i],latinds,loninds)