import settings_interface as si 
import numpy as np
import re

def HEMCOsetup(foldername=None,returnStartEndDict=False):
	spc_config = si.getSpeciesConfig()
	if foldername:
		hemco_config_path = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}/ensemble_runs/{foldername}/"
	else:
		hemco_config_path = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}/template_run/"
	with open(hemco_config_path+'HEMCO_Config.rc') as f:
		lines = f.readlines()
	linenums = np.arange(0,len(lines))
	switchnames,switchcategories,switchboolvals = getHEMCOSwitches(spc_config, lines)
	startlinedict,endlinedict,scaleFactorLineAdd = getHEMCOLines(lines,switchnames,switchcategories,switchboolvals)
	speciesloc = getValidSpeciesAddresses(lines,linenums,startlinedict,endlinedict)
	to_return = [spc_config,lines,hemco_config_path,scaleFactorLineAdd,speciesloc]
	if returnStartEndDict:
		to_return = to_return + [startlinedict,endlinedict,linenums]
	return to_return

def getHEMCOSwitches(spc_config, lines):
	switchnames = []
	switchcategories = []
	switchboolvals = []
	emisregion = 0 #0 for first part (master inputs), 1 for emissions, 2 for non-emissions, 3 for extensions
	inGFED = False #If inGFED, we want to check for the arrow switches; not implemented
	inFINN = False #Similarly if in FINN we want to check for arrow switches; not implemented
	for line in lines:
		#Leave after end of opening switches
		if line.startswith('### END SECTION EXTENSION SWITCHES'):
			break
		elif line.startswith('# ----- REGIONAL'):
			emisregion = 1
		elif line.startswith('# ----- NON-EMISSIONS'):
			emisregion = 2
		elif line.startswith('100'): #Start of HEMCO emissions plugins, which are on off not true false
			emisregion = 3
		if emisregion<3:
			if line.startswith('    --> '): #Base emission component switch
				switchname = line[8:31].strip()
				if line[39:44].strip()=="true":
					boolval = True
				elif line[39:44].strip()=="false":
					boolval = False
				else:
					boolval = None
				switchnames.append(switchname)
				switchcategories.append(emisregion)
				switchboolvals.append(boolval)
		else:
			if line[0:3].isnumeric():  #This row will have on or off.
				switchname = line[8:31].strip()
				if line[33:36].strip()=="on":
					boolval = True
				elif line[33:36].strip()=="off":
					boolval = False
				else:
					boolval = None
				switchnames.append(switchname)
				switchcategories.append(emisregion)
				switchboolvals.append(boolval)
	if spc_config['sim_name']=='fullchem':
		paranoxind = switchnames.index('ParaNOx')
		#Paranox also has a not.Paranox we have to worry about
		paranoxbool = switchboolvals[paranoxind]
		switchnames.append('.not.ParaNOx')
		switchcategories.append(3)
		switchboolvals.append(not paranoxbool)
	return [switchnames,switchcategories,switchboolvals]

def getHEMCOLines(lines,switchnames,switchcategories,switchboolvals):
	switchopenparen = [f'((({name}' for name,boolean,category in zip(switchnames,switchboolvals,switchcategories) if boolean and ((category==1) or (category==3) or (name=='GC_BCs'))]
	switchcloseparen = [f'))){name}' for name,boolean,category in zip(switchnames,switchboolvals,switchcategories) if boolean and ((category==1) or (category==3) or (name=='GC_BCs'))]
	names_on = [name for name,boolean,category in zip(switchnames,switchboolvals,switchcategories) if boolean and ((category==1) or (category==3) or (name=='GC_BCs'))]
	startlinedict = {name:None for name in names_on}
	endlinedict = {name:None for name in names_on}
	#Determine which parts of ship we are skipping
	ship_on = 'SHIP' in names_on
	paranox_on = 'ParaNOx' in names_on
	#Ship has other emissions under it; need to track if we are in that environment
	in_ship = False
	in_ParaNOx = False
	in_notParaNOx = False
	#Scaling factors split by emissions; we need to track if in that environment
	in_sf = False
	#We want to get the first line within Emissions tag under scaling factors;
	#we will want to add our additional scaling factor details here.
	scaleFactorLineAdd = -1 #Line to add
	found_LineAdd = False #Toggle when we find
	for num,line in enumerate(lines):
		if line.startswith('### BEGIN SECTION MASKS'):
			break #Nothing else we need
		if line.startswith('### BEGIN SECTION SCALE FACTORS'):
			in_sf = True #In scale factors now; this takes us all the way to Masks, where we break
		if '(((SHIP\n' in line:
			in_ship = True
		if '(((ParaNOx\n' in line:
			in_ParaNOx = True
		if '(((.not.ParaNOx\n' in line:
			in_notParaNOx = True
		if ')))SHIP\n' in line:
			in_ship = False
		if ')))ParaNOx\n' in line:
			in_ParaNOx = False
		if '))).not.ParaNOx\n' in line:
			in_notParaNOx = False
		if in_ship and (not ship_on):
			continue #Skip line if we are in SHIP and SHIP is off
		if in_ParaNOx and (not paranox_on):
			continue #Skip line if we are in ParaNOx and ParaNOx is off
		if in_notParaNOx and paranox_on:
			continue #Skip line if we are in .not.ParaNOx and ParaNOx is on
		#First emssions tag in scaling factor section is where we will add settings
		if in_sf and ('(((EMISSIONS\n' in line) and (not found_LineAdd):
			scaleFactorLineAdd = num+1 #This is where we insert our stuff.
			found_LineAdd = True
		#Toggle on and off what area of SHIP we are in
		for op in switchopenparen:
			if f'{op}\n' in line:
				name = op[3::]
				if in_sf:
					name = name+'+SCALING_FACTORS'
				if in_ParaNOx and (name != 'ParaNOx'):
					name = name+'+ParaNOx'
				if in_notParaNOx and (name != '.not.ParaNOx'):
					name = name+'+notParaNOx'
				startlinedict[name]=num
		for op in switchcloseparen:
			if f'{op}\n' in line:
				name = op[3::]
				if in_sf:
					name = name+'+SCALING_FACTORS'
				if in_ParaNOx and (name != 'ParaNOx'):
					name = name+'+ParaNOx'
				if in_notParaNOx and (name != '.not.ParaNOx'):
					name = name+'+notParaNOx'
				endlinedict[name]=num
	#Remove entries we did not find; right now just GFED/FINN extension, Volcano, and inorganic iodine. Seems fine.
	startlinedict = {key:val for key, val in startlinedict.items() if val}
	endlinedict = {key:val for key, val in endlinedict.items() if val}
	return [startlinedict,endlinedict,scaleFactorLineAdd]

#Get the line numbers and corresponding section names 
def getValidSpeciesAddresses(lines,linenums,startlinedict,endlinedict):
	#Get names without scaling factors
	section_names = list(startlinedict.keys())
	section_names = [name for name in section_names if not ('+SCALING_FACTORS' in name)]
	#Dictionary addressing line numbers where species are found
	speciesloc = {}
	for name in section_names:
		start = startlinedict[name]
		end = endlinedict[name]
		lines_to_read = lines[start:end]
		linenums_to_read = linenums[start:end]
		for line,num in zip(lines_to_read,linenums_to_read):
			splitline = line.split()
			if len(splitline)==12: #This is a line including emissions settings we might want to modify
				if not splitline[0].isnumeric():
					continue #Usually this means it is commented out.
				species = splitline[8]
				#If we already have found this species, add line number to the list in dictionary
				if species in list(speciesloc.keys()):
					speciesloc[species].append(num)
				else: #Otherwise add a new entry to dictionary
					speciesloc[species] = [num]
	#Remove duplicates
	for species in list(speciesloc.keys()):
		speciesloc[species] = list(set(speciesloc[species]))
	return speciesloc

#Add ScalIDs to lines that we would like to scale. Doesn't differentiate if multiple copies of the same species are present.
#The user is expected to modify
def addScalingFactorNumbers(spc_config,speciesloc,lines):
	species_to_add = spc_config["CONTROL_VECTOR_EMIS"].values()
	specnames = spc_config["CONTROL_VECTOR_EMIS"].keys()
	species_scalid = {} #Dictionary with species emission and scale factor id
	specieskeyval = 700
	for species,specname in zip(species_to_add,specnames):
		print(f'Appending scaling factor IDs for {specname} in HEMCO_Config.')
		if type(species) is list: #Combine these species into one list of line nums to change
			linenums_to_modify = []
			for subspec in species:
				if subspec in speciesloc.keys():
					linenums_to_modify = linenums_to_modify + speciesloc[subspec]
				else:
					print(f'Warning: did not detect species {subspec} in active HEMCO inventories.')
			linenums_to_modify = list(set(linenums_to_modify)) #Remove duplicates
		else:
			if species in speciesloc.keys():
				linenums_to_modify = speciesloc[species]
			else:
				linenums_to_modify = []
				print(f'Warning: did not detect species {species} in active HEMCO inventories.')
		for num in linenums_to_modify: #Go through these line nums and add new scale id
			line = lines[num]
			scalid = line.split()[9]
			notwhitespaceloc = [m.start() for m in re.finditer(r'\S+', line)]
			pre_scalid = line[notwhitespaceloc[0]:notwhitespaceloc[9]]
			line_scalid_section = line[notwhitespaceloc[9]:notwhitespaceloc[10]]
			post_scalid = line[notwhitespaceloc[10]::]
			if scalid == '-':
				#If scalid is - then replace with key
				newline_scalid_section = line_scalid_section.replace(f'{scalid}',f'{specieskeyval}')
			else:
				#Replace scalid with new set which includes addtional key
				newline_scalid_section = line_scalid_section.replace(f'{scalid}',f'{scalid}/{specieskeyval}')
			newline = pre_scalid+newline_scalid_section+post_scalid
			lines[num] = newline #Overwrite line.
		species_scalid[specname] = specieskeyval #Add scalid to dictionary
		specieskeyval+=1 #Increment
		return [lines,species_scalid]

def addScalingFactorFile(spc_config,lines,species_scalid,hemco_config_path,scaleFactorLineAdd):
	#Get range of years from settings.
	if spc_config["DO_ENS_SPINUP"]=="true":
		yearrange = f'{int(np.floor(int(spc_config["ENS_SPINUP_START"])/10000))}-{int(np.floor(int(spc_config["END_DATE"])/10000))}'
	else:
		yearrange = f'{int(np.floor(int(spc_config["START_DATE"])/10000))}-{int(np.floor(int(spc_config["END_DATE"])/10000))}'
	for species in list(species_scalid.keys()):
		print(f'Linking scaling factor files for {species} in HEMCO_Config.')
		scalid = species_scalid[species]
		#Exact and forced timestand. So need the range to include the timestamp on the netCDF
		scalefactor_line_to_add = f'{scalid} ASSIM_{species}  {hemco_config_path}{species}_SCALEFACTOR.nc Scalar {yearrange}/1-12/1-31/0-23 RF xy 1  1\n'
		lines.insert(scaleFactorLineAdd,scalefactor_line_to_add) #insert line at appropriate location
		scaleFactorLineAdd += 1
	return lines

#Update boundary condition section if applicable.
def updateBoundaryConds(spc_config,lines,linenums,startlinedict,endlinedict):
	if spc_config['NEST']=='T':
		print('Updating boundary conditons in HEMCO_Config.')
		lines_to_delete = linenums[(startlinedict['GC_BCs']+1):(endlinedict['GC_BCs'])]
		for i in lines_to_delete:
			lines[i]='' #delete lines
		lines[startlinedict['GC_BCs']+1]=f"* BC_  {spc_config['BC_FILES']} SpeciesBC_?ADV?  1980-2021/1-12/1-31/* EFY xyz 1 * - 1 1"
	else:
		print('Skipping boundary condition update in HEMCO_Config.')
	return lines

def ensureRestartsInTLD(lines,linenums,startlinedict,endlinedict):
	restart_lines = linenums[(startlinedict['GC_RESTART']+1):(endlinedict['GC_RESTART'])]
	for i in restart_lines:
		lines[i] = lines[i].replace("./Restarts/","./")
	return lines

def updateOHforCH4(spc_config,lines):
	if (spc_config["Extensions"]["CH4"] == "True") and (spc_config["USE_CUSTOM_CH4_OH_ENTRY"] == "True"):
		ohval = spc_config["CUSTOM_CH4_OH_ENTRY"]
		for num,line in enumerate(lines):
			if line.startswith("* GLOBAL_OH"):
				lines[num] = f"{ohval}\n"
				break
	return lines

#Write HEMCO Config file. Save different names if it's spinup or nature (only update background) or main (adds scaling factors).
def writeHEMCOConfig(hemco_config_path,lines,spinup_or_nature = False):
	if spinup_or_nature:
		with open(hemco_config_path+'HEMCO_Config_SPINUP_NATURE_TEMPLATE.rc', 'w') as f:
			for line in lines:
				f.write(line)
	else:
		with open(hemco_config_path+'HEMCO_Config.rc', 'w') as f:
			for line in lines:
				f.write(line)


def fullWorkflow(useString, foldername=None):
	spc_config,lines,hemco_config_path,scaleFactorLineAdd,speciesloc,startlinedict,endlinedict,linenums = HEMCOsetup(foldername,returnStartEndDict=True)
	lines = updateOHforCH4(spc_config,lines) #ALL use this update OH script, which only applies if applicable
	lines = ensureRestartsInTLD(lines,linenums,startlinedict,endlinedict) #make sure we are reading restart from TLD, not Restarts folder.
	#self.updateBoundaryConds() #Currently BCs updated by setup Ensemble; no need.
	if useString == 'ENSEMBLE':
		#Add updates for scaling factors
		lines,species_scalid=addScalingFactorNumbers(spc_config,speciesloc,lines)
		lines = addScalingFactorFile(spc_config,lines,species_scalid,hemco_config_path,scaleFactorLineAdd)
		writeHEMCOConfig(hemco_config_path,lines,False)
	elif useString == 'NATURE':
		writeHEMCOConfig(hemco_config_path,lines,True)


