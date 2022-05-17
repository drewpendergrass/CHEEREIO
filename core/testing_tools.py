import letkf_utils as lu
import numpy as np

assimilator = lu.Assimilator('20190108_0000',2,1) #we aren't using ens or core num here.

def walkThroughAssimilation(assim,latind=65,lonind=24): #default is a point in northern California, arbitrary 
	print('*************************************************')
	print('*******ASSIMILATION WALKTHROUGH, VERBOSE*********')
	print('*************************************************')
	assim.verbose=2
	print(f'Coordinates for this walkthrough: {(assim.gt[1].getLat()[latind],assim.gt[1].getLon()[lonind])}')
	assim.prepareMeansAndPerts(latind,lonind)
	print(f'ybar has value {assim.ybar_background}')
	print(f'ydiff has value {assim.ydiff}')
	assim.makeR(latind,lonind)
	sqrtdiag = np.sqrt(np.diag(assim.R))
	print(f'The square root of the diagonal of R is {sqrtdiag}')
	print(f'Relative to ybar, that diagonal is value {sqrtdiag/assim.ybar_background}')
	assim.makeC()
	assim.makePtildeAnalysis()
	assim.makeWAnalysis()
	assim.makeWbarAnalysis()
	assim.adjWAnalysis()
	assim.makeAnalysisCombinedEnsemble()
	analysisSubset,backgroundSubset,analysisPertSubset,backgroundPertSubset = assim.getAnalysisAndBackgroundColumn(latind,lonind,doBackground=True,doPerts=True)
	print(f'Background column has dimension {np.shape(backgroundSubset)} and value {backgroundSubset}')
	print(f'Analysis column has dimension {np.shape(analysisSubset)} and value {analysisSubset}')
	meananalysis = np.mean(analysisSubset,axis=1)
	meanbackground = np.mean(backgroundSubset,axis=1)
	print(f'Ensemble mean analysis minus ens mean background has value {meananalysis-meanbackground}')
	print(f'This represents a percent difference of {100*((meananalysis-meanbackground)/meanbackground)}%')
	dofs = assim.calculateDOFS(analysisPertSubset,backgroundPertSubset)
	print(f'DOFS has value {dofs}')
	analysisSubsetAdjusted = assim.applyAnalysisCorrections(analysisSubset,backgroundSubset,latind,lonind)
	meananalysisSubsetAdjusted = np.mean(analysisSubsetAdjusted,axis=1)
	print(f'Old analysis ensemble mean for column was {meananalysis}.')
	print(f'After postprocessing, new analysis ensemble mean for column is {meananalysisSubsetAdjusted}.')
	print(f'This represents a percent difference of {100*((meananalysisSubsetAdjusted-meananalysis)/meananalysis)}%')
