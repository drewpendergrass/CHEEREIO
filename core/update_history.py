import settings_interface as si 
import sys
import numpy as np

class HISTORY_Translator():
	#Constructor. If foldername is none, open template HISTORY.rc. Otherwise open ensemble member
	def __init__(self,foldername=None):
		self.spc_config = si.getSpeciesConfig()
		if not foldername:
			#Default to template HEMCO config
			self.historyrc_path = f"{self.spc_config['MY_PATH']}/{self.spc_config['RUN_NAME']}/template_run/"
		else:
			#Otherwise go to Control, or specific ensemble member
			if foldername=='CONTROL':
				self.historyrc_path = f"{self.spc_config['MY_PATH']}/{self.spc_config['RUN_NAME']}/control_run/"
			else:
				self.historyrc_path = f"{self.spc_config['MY_PATH']}/{self.spc_config['RUN_NAME']}/ensemble_runs/{foldername}/"
		with open(self.historyrc_path+'HISTORY.rc') as f:
			self.lines = f.readlines()
		self.linenums = np.arange(0,len(self.lines))
	def calcDurFreq(self,isFirst):
		if isFirst=="First":
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
		elif isFirst=="Spinup":
			timestr="'End'"
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
				self.lines[num] = f'  Restart.frequency:          {timestr}\n'
			if line.startswith('  Restart.duration'):
				self.lines[num] = f'  Restart.duration:           {timestr}\n'
				break
	def updateHistoryCollectionsDurationFrequency(self,isSpinup):
		if isSpinup:
			freqstr = self.spc_config['SPINUP_HISTORY_freq']
			durstr = self.spc_config['SPINUP_HISTORY_dur']
		else:
			freqstr = self.spc_config['HISTORY_freq']
			durstr = self.spc_config['HISTORY_dur']
		for num,line in enumerate(self.lines):
			for collection in self.spc_config["HISTORY_collections_to_customize"]:
				if line.startswith(f'  {collection}.frequency'):
					whitespace = " "*(17-len(collection))
					self.lines[num] = f'  {collection}.frequency:{whitespace}{freqstr}\n'
				if line.startswith(f'  {collection}.duration'):
					whitespace = " "*(18-len(collection))
					self.lines[num] = f'  {collection}.duration:{whitespace}{durstr}\n'
	def findSectionLines(self,sectionname):
		counting = False
		startstop = []
		for num,line in enumerate(self.lines):
			if line.startswith(f'  {sectionname}.fields:'):
				counting = True
				startstop.append(num)
				continue
			if counting and line.startswith('::'):
				counting = False
				startstop.append(num)
				break
		return startstop
	def prepLevelEdgeDiags(self):
		if self.spc_config['SaveLevelEdgeDiags']=='True':
			print('Turning on the LevelEdgeDiags collection in HISTORY.rc!')
			for num,line in enumerate(self.lines):
				if "#'LevelEdgeDiags'," in line:
					self.lines[num] = self.lines[num].replace('#','')
					break
		else:
			print('Turning off the LevelEdgeDiags collection in HISTORY.rc!')
	def prepStateMet(self):
		if self.spc_config['SaveStateMet']=='True':
			print('Turning on the StateMet collection in HISTORY.rc!')
			for num,line in enumerate(self.lines):
				if "#'StateMet'," in line:
					self.lines[num] = self.lines[num].replace('#','')
					break
		else:
			print('Turning off the StateMet collection in HISTORY.rc!')
	def makeSectionString(self,species,isFirst,sectionname):
		if sectionname=='SpeciesConc':
			if isFirst:
				startstring = "  SpeciesConc.fields:         "
			else:
				startstring = "                              "
			secondstring = "'SpeciesConc_"
			endwhitespacecount=18-len(species)
		elif sectionname=='LevelEdgeDiags':
			if isFirst:
				startstring = "  LevelEdgeDiags.fields:      "
			else:
				startstring = "                              "
			secondstring = "'"
			endwhitespacecount=30-len(species)
		elif sectionname=='StateMet':
			if isFirst:
				startstring = "  StateMet.fields:            "
			else:
				startstring = "                              "
			secondstring = "'"
			endwhitespacecount=30-len(species)
		finalstring = "',\n"
		return startstring+secondstring+species+(' '*endwhitespacecount)+finalstring
	def customizeSection(self,sectionname):
		print(f'Overwriting default {sectionname} collection in HISTORY.rc.')
		startstop = self.findSectionLines(sectionname)
		del self.lines[startstop[0]:startstop[1]]
		isFirst = True
		for count, species in enumerate(self.spc_config[f'History{sectionname}ToSave']):
			print(f'Adding {species} to the {sectionname} collection in HISTORY.rc.')
			self.lines.insert(startstop[0]+count, self.makeSectionString(species,isFirst,sectionname))
			isFirst = False
	def customizeAllSections(self):
		self.customizeSection('SpeciesConc')
		if self.spc_config['SaveLevelEdgeDiags']=='True':
			self.customizeSection('LevelEdgeDiags')
		if self.spc_config['SaveStateMet']=='True':
			self.customizeSection('StateMet')
	def writeHistoryConfig(self):
		with open(self.historyrc_path+'HISTORY.rc', 'w') as f:
			for line in self.lines:
				f.write(line)
		print('HISTORY.rc saved successfully. Please double check that all the collections match what you want to save!')

settingsstr = str(sys.argv[1])

if settingsstr=="SETCONTROL":
	trans = HISTORY_Translator(foldername='CONTROL')
else:
	trans = HISTORY_Translator()

if settingsstr=="TEMPLATEDIR":
	trans.prepLevelEdgeDiags()
	trans.prepStateMet()
	trans.updateRestartDurationFrequency(isFirst="First")
	trans.customizeAllSections()
	trans.updateHistoryCollectionsDurationFrequency(isSpinup=False)
elif settingsstr=="SPINUP":
	trans.prepLevelEdgeDiags()
	trans.prepStateMet()
	trans.updateRestartDurationFrequency(isFirst="Spinup")
	trans.customizeAllSections()
	trans.updateHistoryCollectionsDurationFrequency(isSpinup=True)
elif settingsstr=="PREPMAIN":
	trans.updateHistoryCollectionsDurationFrequency(isSpinup=False)
elif settingsstr=="UPDATEDURFREQ":
	trans.updateRestartDurationFrequency(isFirst="Midrun")
elif settingsstr=="SETCONTROL":
	trans.updateRestartDurationFrequency(isFirst="Midrun")
	trans.updateHistoryCollectionsDurationFrequency(isSpinup=False)


trans.writeHistoryConfig()