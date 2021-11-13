import toolbox as tx 

class HISTORY_Translator():
	#Constructor. If foldername is none, open template HISTORY.rc. Otherwise open ensemble member
	def __init__(self,foldername=None,testing=False):
		self.spc_config = tx.getSpeciesConfig(testing)
		if foldername:
			self.foldername = foldername
		else:
			self.foldername = ''
		if not foldername:
			#Default to template HEMCO config
			self.historyrc_path = f"{self.spc_config['MY_PATH']}/{self.spc_config['RUN_NAME']}/template_run/"
		else:
			#Otherwise go to specific ensemble member
			self.historyrc_path = f"{self.spc_config['MY_PATH']}/{self.spc_config['RUN_NAME']}/ensemble_runs/{foldername}/"
		with open(self.historyrc_path+'HISTORY.rc') as f:
			self.lines = f.readlines()
		self.linenums = np.arange(0,len(self.lines))
	def findSpecConcLines(self):
		counting = False
		startstop = []
		for num,line in enumerate(self.lines):
			if line.startswith('  SpeciesConc.fields:'):
				counting = True
				startstop.append(num)
				continue
			if counting and line.startswith('::'):
				counting = False
				startstop.append(num)
				break
		self.SpecConcStartStopLines = startstop
	def prepLevelEdgeDiags(self):
		if self.spc_config['SaveLevelEdgeDiags']=='True':
			print('Turning on the LevelEdgeDiags collection in HISTORY.rc!')
			for num,line in enumerate(self.lines):
				if "#'LevelEdgeDiags'," in line:
					self.lines[num] = self.lines[num].replace('#','')
					break
		else:
			print('Turning off the LevelEdgeDiags collection in HISTORY.rc!')
	def makeSpecConcString(self,species,isFirst):
		if isFirst:
			startstring = "  SpeciesConc.fields:         "
		else:
			startstring = "                              "
		secondstring = "SpeciesConc_"
		endwhitespacecount=18-len(species)
		finalstring = "',\n"
		return startstring+secondstring+species+(' '*endwhitespacecount)+finalstring
	def customizeSpecConc(self):
		self.findSpecConcLines()
		print('Overwriting default SpeciesConc collection in HISTORY.rc.')
		del self.lines[self.SpecConcStartStopLines[0]:self.SpecConcStartStopLines[1]]
		isFirst = True
		for species in self.spc_config['HistorySpecConcToSave']:
			print(f'Adding {species} to the SpeciesConc collection in HISTORY.rc.')
			self.lines.insert(self.SpecConcStartStopLines[0], self.makeSpecConcString(species,isFirst))
			isFirst = False
	def writeHistoryConfig(self):
		with open(self.historyrc_path+'HISTORY.rc', 'w') as f:
			for line in self.lines:
				f.write(line)
		print('HISTORY.rc saved successfully. Please double check that all the collections match what you want to save!')

trans = HISTORY_Translator()
trans.prepLevelEdgeDiags()
trans.customizeSpecConc()
trans.writeHistoryConfig()