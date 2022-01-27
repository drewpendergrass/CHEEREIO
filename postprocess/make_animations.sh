#!/bin/bash

eval "$(conda shell.bash hook)"

# Name for this run
RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"
animation_fps_scalingfactor="$(jq -r ".animation_fps_scalingfactor" ../ens_config.json)"
animation_fps_concentrations="$(jq -r ".animation_fps_concentrations" ../ens_config.json)"
# Path where you want to set up assimilation code and run directories
MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
pp_dir="${MY_PATH}/${RUN_NAME}/postprocess" 
concpp="${pp_dir}/controlvar_pp.nc"


controlvec=( $(jq -r ".CONTROL_VECTOR_CONC[]" ../ens_config.json) )
emisvec=( $(jq -r '.CONTROL_VECTOR_EMIS | keys[] as $k | "\($k)"' ../ens_config.json) )

conda activate $(jq -r ".AnimationEnv" ../ens_config.json) 

## now loop through the control vector array
for i in "${controlvec[@]}"
do
	python animator.py "${concpp}" "SpeciesConc_${i}" "all" "${pp_dir}/SpeciesConc_${i}" "${animation_fps_concentrations}"
done

for i in "${emisvec[@]}"
do
	python animator.py "${pp_dir}/${i}_SCALEFACTOR.nc" "Scalar" "all" "${pp_dir}/SCALEFACTOR_${i}" "${animation_fps_scalingfactor}"
done

conda deactivate

exit 0