# ITCZstate_JGRA
JGRA ITCZ states paper scripts README

# Necessary procedures before running scripts to produce figures

Bilinearly daily precipitation data to a common 1 x 1 degree regular grid using Climate Data Operators (CDO, https://code.mpimet.mpg.de/projects/cdo) *or a programming language of choice(e.g., “cdo remapbil,r360x181 input_file.nc output_file.nc”)

Shell script obs and reanalyses: “ITCZ_states_absthresh_obs_reana_savenc.sh”

Shell script for CMIP6: “ITCZ_states_absthresh_CMIP6_savenc.sh”

Python file to calculate ITCZ states: “ITCZ_states_absthresh_savenc.py”

# Python scripts that produce figures in the main paper

Fig. 1 – “Annual_precip_bias.ipynb” – produces maps of the global annual mean precipitation in IMERG observations and CMIP6 models

Fig. 3 – “ITCZst_examples.ipynb” – produces example days of each ITCZ state in IMERG observations

Figs. 4, 7, 8 and 9 – “ITCZ_states_monthly_mean_stats.ipynb” – produces mean ITCZ states vs. month in CMIP6 models

Fig. 5 – “ITCZst_comparison_maps_obs_reana.ipynb” – produces maps of mean daily precipitation for each ITCZ state in observations and reanalyses; also computes statistics in Table 2.

Fig. 6 – “Interann_ITCZ_state_obs_reana_CMIP6_savenc.ipynb” – calculates monthly frequencies of ITCZ state in observations, reanalyses, and CMIP6 models, then plots standard deviation of ITCZ in month for observations and reanalyses

Fig. 10 – “ITCZst_comparison_maps_obs_CMIP6.ipynb” – produces maps of mean daily precipitation for each ITCZ state in observations and CMIP6 model groups; also computes statistics in Table 2.

Fig. 11 – “Interann_ITCZ_state_variability_CMIP6_groups.ipynb” – produces plots of standard deviation of ITCZ in month for observations, reanalyses, and CMIP6 models

Fig. 12 – “Precip_thresh_obs_reana_CMIP6.ipynb” – compute statistical distributions of 90th percentile precipitation thresholds for dITCZs, sITCZs, and nITCZs in observations, reanalyses, and CMIP6 models

# Python scripts that produce figures in the supporting information file

Fig. S1 – “ITCZ_state_errors_allbasins_supporting.ipynb” – plots ITCZ state error averaged over all ITCZ states and months in reanalyses and CMIP6 models for the east Pacific, central Pacific, Atlantic, and Indian Ocean regions

Fig. S2 – “Obs_reanalysis_event_duration_supporting.ipynb” – produces interquartile range of the number of days of each ITCZ state in observations and reanalyses

Fig. S3 – “Obs_ITCZ_states_abs_thresh_supporting.ipynb” – ITCZ state as a function of absolute precipitation threshold, from 1 to 9.5 mm/day in IMERG observations
