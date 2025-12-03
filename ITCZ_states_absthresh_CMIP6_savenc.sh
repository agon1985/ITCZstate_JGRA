#!/bin/bash
# This is a shell script to run the ITCZ states algorithm (which is a Python script) for multiple datasets at once, 
# given they have the same file structure of the precipitation field. In this case, we are running 25 CMIP6 models.

#Author: Alex Gonzalez
#Date: November 2025
#Paper: JGR-Atmospheres

# Navigate to the directory where the Python script is located (optional)
cd /home/alex.gonzalez/CMIP6_ITCZ_state_project/scripts/methods_paper/revision/final/

acron_array=('BCC-CSM2-MR'
            'CAMS-CSM1-0'
            'CESM2'
            'CESM2-WACCM'
            'CMCC-CM2-HR4'
            'CMCC-CM2-SR5'
    	    'CMCC-ESM2'
            'CNRM-CM6-1-HR'
            'E3SM-1-0'
            'E3SM-2-0'
            'E3SM-2-0-NARRM'
            'EC-Earth3'
            'EC-Earth3-AerChem'
            'EC-Earth3-CC'
            'EC-Earth3-Veg'
            'EC-Earth3-Veg-LR'
            'FGOALS-f3-L'
            'GFDL-CM4'
            'GFDL-ESM4'
            'HadGEM3-GC31-MM'
            'MPI-ESM1-2-HR'
            'MRI-ESM2-0'
            'NorESM2-MM'
            'SAM0-UNICON'
    	    'TaiESM1')

main_folder='/srv/lss/CMIP6_precip_1deg/' # main folder
nc_name='pr'	 # netcdf name of precipitation field
unit_const=86400 # to convert to mm day^{-1}

# original netcdf names for latitude and longitude
tim_name='time'
lat_name='lat'
lon_name='lon'

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
#region='WestPac'
#region='CentPac'

region_str='East_Pacific_Ocean_(90-135W)'
#region_str='Atlantic_Ocean_(0-45W)'
#region_str='Indian_Ocean_(50-95E)'
#region_str='West_Pacific_Ocean_(165E-150W)'
#region_str='Central_Pacific_Ocean_(135W-180W)'

# East Pac
lon0=-135
lon1=-90

# Atlantic
#lon0=-45
#lon1=0

# Indian Ocean
#lon0=50
#lon1=95

# West Pac
#lon0=165
#lon1=210

# Central Pac
#lon0=-180
#lon1=-135
#----------------------------------------------------------------#
#----------------------------------------------------------------#

ndataset=${#acron_array[@]}
echo "There are $ndataset total datasets."

for item in $(seq 0 $(($ndataset - 1))); do
    echo "This is the dataset called ${acron_array[$item]}."
    folder_name=$main_folder${acron_array[$item]}'/'
    file_name='pr_'${acron_array[$item]}'_1980-2014_1deg_180.nc'
    #echo $folder_name
    #echo $file_name
    python ITCZ_state_algorithm_absthresh_savenc.py $abs_thresh ${acron_array[$item]} $folder_name $file_name $nc_name $tim_name $lat_name $lon_name $unit_const $region $region_str $lon0 $lon1

done	    
