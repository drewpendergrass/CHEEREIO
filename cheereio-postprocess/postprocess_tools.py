import sys
import os
import json

with open('../ens_config.json') as f:
	ens_config_data = json.load(f)

sys.path.append(os.path.abspath(f'{ens_config_data['ASSIM_PATH']}/core/'))

import letkf_utils as lu