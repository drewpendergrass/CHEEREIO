import settings_interface as si 
import ruamel.yaml
import sys
import os.path

spc_config = si.getSpeciesConfig()
obspack_code = spc_config["ACTIVATE_OBSPACK"] == 'true'
obspack_file = spc_config["obspack_gc_input_file"]
obspack_path = spc_config["gc_obspack_path"]
path = str(sys.argv[1])

yaml = ruamel.yaml.YAML()
# yaml.preserve_quotes = True
with open(path) as fp:
	data = yaml.load(fp)

#switch obspack to user setting
data['extra_diagnostics']['obspack']['activate'] = obspack_code
data['extra_diagnostics']['obspack']['input_file'] = f'{obspack_path}/{obspack_file}'

with open(f'{path}/geoschem_config.yml', "w") as f:
	yaml.dump(data, f)