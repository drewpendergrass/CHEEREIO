Tips for nested grid simulations
==========

CHEEREIO supports inversions over nested grid regions (e.g. North America or China, as opposed to the world). This feature is not yet seamlessly supported in CHEEREIO, but neverless users have gotten nested simulation inversions to work and have gotten interesting results. Here are some community-sourced tips to get nested grid simulations to work:

Users doing nested grid inversions are encouraged to ensure consistency across three specific parts of the CHEEREIO setup, especially for AS-nested setups:

#. The latitude/longitude range in the restart file.
#. The AS range defined in core/toolbox.py (i.e. the makeLatLonGridWithMask() function)
#. The filtering of observations by latitude/longitude in the observation operator, specifically in the apply_filters function. In particular, observations outside of the nested domain should be discarded.

Per Huiru Zhong (PKU): To simulate nested assimilation for eastern China, I made the following adjustments in the configuration:

#. Defined AS in the ens_config file for the nested setup
#. Generated a restart file for eastern China with grid settings lon-range[98,135] and lat-range[18,52]
#. Modified the AS region’s longitude and latitude in toolbox.py for AS_MERRA2 (lon = np.arange(98.125,135,0.625) and lat = np.arange(18,52, 0.5))
#. Added filters in CHEEREIO/extensions/TROPOMI_CH4_extension.json and included additional lines of code in observation_operator.py to filter observations by longitude and latitude

Per Huiru Zhong (PKU): When these areas are carefully matched, emissions should not be reduced to zero. However, I overlooked some of these steps recently while submitting new tasks when modifying boundary conditions, which caused the emissions to zero out. I suspect that this issue might be due to a mismatch between the dimensions updated with each assimilation and the actual dimensions read by HEMCO, leading to significant discrepancies between simulations and observations.

Thanks to Yifan Li (Tsingua Shenzhen), Huiru Zhong (Peking University), and Yunxiao Tang (Harvard) for their help in compiling these tips. I plan to enforce these settings in later versions of CHEEREIO, streamlining nested grid inversions.