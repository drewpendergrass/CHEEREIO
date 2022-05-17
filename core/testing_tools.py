import letkf_utils as lu
import numpy as np
import json

#Makes an assimilator object
def makeAssimilator():
	assimilator = lu.Assimilator('20190108_0000',2,1) #we aren't using ens or core num here.
	return assimilator

#Overrides settings without modifying ens_config. If overwrite is true, it deletes previous adjustments
def overrideSettings(settings_to_override, overwrite = False)
	with open('../settings_to_override.json') as f:
		over_data = json.load(f)
	if overwrite:
		over_data = settings_to_override #Flat overwrite
	else: #If not overwriting, loop through and add keys individually
		for key in list(settings_to_override.keys()):
			if key != 'override':
				over_data[key] = settings_to_override[key]
	over_data['override'] = "True" #Activate overriding capability
	with open('../settings_to_override.json', 'w') as f:
		json.dump(over_data, f, ensure_ascii=False, indent=4)

#Walks through with extensive print statements an assimilation cycle
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
	meananalysis = np.mean(analysisSubset,axis=1)
	meanbackground = np.mean(backgroundSubset,axis=1)
	print(f'Background column has dimension {np.shape(backgroundSubset)} and ensemble mean value {meanbackground}')
	print(f'Analysis column has dimension {np.shape(analysisSubset)} and ensemble mean value {meananalysis}')
	print(f'Ensemble mean analysis minus ens mean background has value {meananalysis-meanbackground}')
	print(f'This represents a percent difference of {100*((meananalysis-meanbackground)/meanbackground)}%')
	dofs = assim.calculateDOFS(analysisPertSubset,backgroundPertSubset)
	print(f'DOFS has value {dofs}')
	analysisSubsetAdjusted = assim.applyAnalysisCorrections(analysisSubset,backgroundSubset,latind,lonind)
	meananalysisSubsetAdjusted = np.mean(analysisSubsetAdjusted,axis=1)
	print(f'Old analysis ensemble mean for column was {meananalysis}.')
	print(f'After postprocessing, new analysis ensemble mean for column is {meananalysisSubsetAdjusted}.')
	print(f'This represents a percent difference of {100*((meananalysisSubsetAdjusted-meananalysis)/meananalysis)}%')


