import toolbox as tx
import sys
import numpy as np

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
	def calcDurFreq(self,isFirst):
		if isFirst:
			ASSIM_START_DATE = self.spc_config['ASSIM_START_DATE']
			assim_start_year = int(ASSIM_START_DATE[0:4])
			assim_start_month = int(ASSIM_START_DATE[4:6])
			assim_start_day = int(ASSIM_START_DATE[6:8])
			START_DATE = self.spc_config['START_DATE']
			start_year = int(START_DATE[0:4])
			start_month = int(START_DATE[4:6])
			start_day = int(START_DATE[6:8])
			monthsapart = (assim_start_year-start_year)*12+(assim_start_month-start_month)
			yeardiff = int(np.floor(monthsapart/12))
			monthdiff = monthsapart%12
			daydiff = assim_start_day-start_day
			if daydiff<0:
				print('This date configuration not supported; defaulting to "End" restart setting')
				timestr="'End'"
			else:
				yearstr = str(yeardiff).zfill(4)
				monthstr = str(monthdiff).zfill(2)
				daystr = str(daydiff).zfill(2)
				timestr = f'{yearstr}{monthstr}{daystr} 000000'
		else:
			ASSIM_TIME = int(self.spc_config['ASSIM_TIME'])
			assim_days = int(np.floor(ASSIM_TIME/24))
			assim_hours = ASSIM_TIME%24
			daystr = str(assim_days).zfill(2)
			hourstr = str(assim_hours).zfill(2)
			timestr = f'000000{daystr} {hourstr}0000'
		return timestr
	def updateRestartDurationFrequency(self, isFirst):
		timestr = self.calcDurFreq(isFirst)
		for num,line in enumerate(self.lines):
			if line.startswith('  Restart.frequency'):
				self.lines[num] = f'  Restart.frequency:          {timestr},\n'
			if line.startswith('  Restart.duration'):
				self.lines[num] = f'  Restart.duration:           {timestr},\n'
				break
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

settingsstr = str(sys.argv[1])
trans = HISTORY_Translator()

if settingsstr=="TEMPLATEDIR":
	trans.prepLevelEdgeDiags()
	trans.updateRestartDurationFrequency(isFirst=True)
	trans.findSpecConcLines()
	trans.customizeSpecConc()
elif settingsstr=="UPDATEDURFREQ":
	trans.updateRestartDurationFrequency(isFirst=False)

trans.writeHistoryConfig()
