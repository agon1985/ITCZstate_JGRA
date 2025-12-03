#!/bin/bash
# This is a shell script to execute a Python file

# Navigate to the directory where the Python script is located (optional)
cd /home/alex.gonzalez/CMIP6_ITCZ_state_project/scripts/methods_paper/revision/final/

acron_array=('IMERG' 'IMERG_mw' 'GPCP' 'ERA5' 'MERRA-2' 'JRA-3Q' 'CFSR')

mainfolder='/srv/lss/'
folder_array=('IMERG_1deg/total_precip/' 'IMERG_1deg/total_precip/' 
	'GPCP_daily/total_precip/' 'ERA5_1deg/total_precip/' 
	'MERRA-2/total_precip/' 'JRA-3Q/total_precip/' 'CFSR/total_precip/')


fname_array=('imerg_1deg_1998-2024_daily.nc' 'imerg_mw_1998-2024_daily.nc' 'gpcp_v01r03_daily_1996-2024_180.nc'
	'era5_tp_1979-2024_1deg_180.nc' 'MERRA2_pr_1980-2025_daily_180_cyclic.nc' 'jra-3q_1979-2024_1deg_daily.nc'
	'cfsr_pr_1979-2025_daily_1deg.nc')
nc_name_array=('precipitation' 'MWprecipitation' 'precip' 'tp' 'PRECTOTCORR' 'tprate1have-sfc-fc-gauss' 'PRATE_L1_Avg_1')

tim_name='time'

lat_name_array=('lat' 'lat' 'latitude' 'latitude' 'lat' 'lat' 'lat')
lon_name_array=('lon' 'lon' 'longitude' 'longitude' 'lon' 'lon' 'lon')
unit_const_array=(1 1 1 24000 86400 86400 86400) 

# mm/day precip threshold
abs_thresh=5.0

#----------------------------------------------------------------#
##### THE USER MUST KNOW THE ORDER OF LATITUDE AND LONGITUDE #####
##### AND ORDER ALL VALUES BELOW CONSISTENT WITH THIS ORDER #####
#----------------------------------------------------------------#

# USER INPUT: 1) region abbreviation, 2) region title name, 
# 3) domain western longitude (lon0), 4) domain eastern longitude (lon1)

region='EastPac'
#region='Atlantic_w'
#region='Indian'
#region='CentPac'
#region='WestPac'

region_str='East_Pacific_Ocean_(90-135W)'
#region_str='Atlantic_Ocean_(0-45W)'
#region_str='Indian_Ocean_(50-95E)'
#region_str='Central_Pacific_Ocean_(135W-180W)'
#region_str='West_Pacific_Ocean_(165E-150W)'

# East Pac
lon0=-135
lon1=-90

# Atlantic
#lon0=-45
#lon1=0

# Indian Ocean
#lon0=50
#lon1=95

# Central Pac
#lon0=-180
#lon1=-135

# West Pac
#lon0=165
#lon1=210
#----------------------------------------------------------------#
#----------------------------------------------------------------#

ndataset=${#acron_array[@]}
echo "There are $ndataset total datasets."

for item in $(seq 0 $(($ndataset - 1))); do
    echo "This is the dataset called ${acron_array[$item]}."
    python ITCZ_state_algorithm_absthresh_savenc.py $abs_thresh ${acron_array[$item]} $mainfolder${folder_array[$item]} ${fname_array[$item]} ${nc_name_array[$item]} $tim_name ${lat_name_array[$item]} ${lon_name_array[$item]} ${unit_const_array[$item]} $region $region_str $lon0 $lon1


done    
