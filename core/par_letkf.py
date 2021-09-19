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
# assimilator.makeR(23,0)
# assimilator.makeC()
# assimilator.makePtildeAnalysis()
# assimilator.makeWAnalysis()
# assimilator.makeWbarAnalysis()
# assimilator.adjWAnalysis()
# assimilator.makeAnalysisCombinedEnsemble()

# natureobs = np.concatenate([assimilator.ObsOp['NO_SATELLITE'].obsinfo.getObsVal(23,0),assimilator.ObsOp['NO2_SATELLITE'].obsinfo.getObsVal(23,0)])
# natvec = assimilator.nature.getStateVector(23,0)

# (assimilator.Ypert_background[0,:]+assimilator.ybar_background[0])/natureobs[0]

# (assimilator.Xpert_background[0,:]+assimilator.xbar_background[0])/natvec[0]
# assimilator.analysisEnsemble[0,:]/natvec[0]

# assimilator.inflation=0.1

# scalar = 0.01
# assimilator.ensemble_numbers = np.array([1,2,3,4,5])
# assimilator.ybar_background = np.array([5])*scalar
# assimilator.Ypert_background = np.array([[-4,-2,-0,2,4]])*scalar
# assimilator.ydiff = np.array([0])*scalar
# assimilator.xbar_background = np.array([5])*scalar*4.43605172e-05
# assimilator.Xpert_background = np.array([[-4,-2,-0,2,4]])*scalar*4.43605172e-05
# assimilator.R=np.array([[1]])*scalar
# assimilator.makeC()
# assimilator.makePtildeAnalysis()
# assimilator.makeWAnalysis()
# assimilator.makeWbarAnalysis()
# assimilator.adjWAnalysis()
# assimilator.makeAnalysisCombinedEnsemble()
# assimilator.analysisEnsemble/(scalar*4.43605172e-05)


# assimilator.ybar_background = np.array([8.38947803484942e-08])
# assimilator.Ypert_background = np.array([[-4.93498701e-08, -4.27698888e-08, -3.61899059e-08, -2.96099234e-08,-2.30299405e-08, -1.64499597e-08, -9.86997612e-09, -3.28999328e-09,3.28999537e-09,  9.86997750e-09,  1.64499614e-08,  2.30299417e-08,2.96099242e-08,  3.61899045e-08,  4.27698883e-08,  4.93498648e-08]])
# assimilator.ydiff = np.array([-2.7964924046621343e-08])
# assimilator.xbar_background = np.array([3.721615761875449e-12])
# assimilator.Xpert_background = np.array([[-2.18918576e-12, -1.89729425e-12, -1.60540296e-12, -1.31351145e-12,-1.02162015e-12, -7.29728642e-13, -4.37837348e-13, -1.45945837e-13, 1.45945891e-13,  4.37837402e-13,  7.29728913e-13,  1.02161999e-12,1.31351150e-12,  1.60540258e-12,  1.89729452e-12,  2.18918560e-12]])
# assimilator.R=np.array([[5.592985630187286e-09]])
# assimilator.makeC()
# assimilator.makePtildeAnalysis()
# assimilator.makeWAnalysis()
# assimilator.makeWbarAnalysis()
# assimilator.adjWAnalysis()
# assimilator.makeAnalysisCombinedEnsemble()
# assimilator.analysisEnsemble/natvec[0]
# (assimilator.Xpert_background+assimilator.xbar_background)/natvec[0]

# scaler = 8.38947803484942e-08
# scaler = 5
# assimilator.ybar_background = np.array([1])*scaler
# assimilator.Ypert_background = np.array([[-0.58823529, -0.50980393, -0.43137256, -0.35294119, -0.27450981,-0.19607846, -0.11764708, -0.03921571,  0.03921573,  0.1176471 ,0.19607848,  0.27450983,  0.3529412 ,  0.43137254,  0.50980392,0.58823522]])*scaler
# #assimilator.ydiff = np.array([0])
# assimilator.ydiff = assimilator.ybar_background / -3
# assimilator.R = np.array([[assimilator.ybar_background[0] / 15]])
# #assimilator.R = np.array([0.05])
# assimilator.xbar_background = assimilator.ybar_background * 4.43605172e-05
# assimilator.Xpert_background = assimilator.Ypert_background * 4.43605172e-05
# assimilator.makeC()
# assimilator.makePtildeAnalysis()
# assimilator.makeWAnalysis()
# assimilator.makeWbarAnalysis()
# assimilator.adjWAnalysis()
# assimilator.makeAnalysisCombinedEnsemble()
# assimilator.analysisEnsemble/(4.43605172e-05*scaler)
# (assimilator.Xpert_background+assimilator.xbar_background)/(4.43605172e-05*scaler)



# assimilator.prepareMeansAndPerts(23,0)
# assimilator.makeR(23,0)

# np.diag(assimilator.NatureHelperInstance.obs_info.errs['NO_SATELLITE'])/assimilator.NatureHelperInstance.obs_info.values['NO_SATELLITE']
# np.diag(assimilator.NatureHelperInstance.obs_info.errs['NO2_SATELLITE'])/assimilator.NatureHelperInstance.obs_info.values['NO2_SATELLITE']
# np.diag(assimilator.NatureHelperInstance.obs_info.errs['NO2_SURFACE'])/assimilator.NatureHelperInstance.obs_info.values['NO2_SURFACE']

# np.diag(assimilator.O.obs_info.errs['NO_SATELLITE'])/assimilator.NatureHelperInstance.obs_info.values['NO_SATELLITE']

# assimilator.prepareMeansAndPerts(23,0)
# assimilator.makeR(23,0)
# assimilator.makeC()
# assimilator.makePtildeAnalysis()
# assimilator.makeWAnalysis()
# assimilator.makeWbarAnalysis()
# assimilator.adjWAnalysis()
# assimilator.makeAnalysisCombinedEnsemble()


# assimilator.prepareMeansAndPerts(23,0)
# assimilator.ybar_background-assimilator.Ypert_background[:,0]
# latinds,loninds = tx.getIndsOfInterest(23,0,testing=True)
# conc4D=assimilator.combineEnsembleForSpecies('NO2')
# obsEns = np.zeros([len(assimilator.ObsOp['NO2_SATELLITE'].obsinfo.getObsVal(23,0)),np.shape(conc4D)[3]])

# for i in range(16):
# 	obsEns[:,i],_,_ = assimilator.ObsOp['NO2_SATELLITE'].H(conc4D[:,:,:,i],latinds,loninds)