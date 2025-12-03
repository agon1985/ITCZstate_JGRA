#!/bin/bash
# This is a shell script to execute a Python file

# Navigate to the directory where the Python script is located (optional)
cd /home/alex.gonzalez/

# Execute the Python script
#python your_script_name.py

my_array=('BCC-CSM2-MR',
            'CAMS-CSM1-0',
            'CESM2',
            'CESM2-WACCM',
            'CMCC-CM2-HR4',
            'CMCC-CM2-SR5',
	    'CMCC-ESM2',
            'CNRM-CM6-1-HR',
            'E3SM-1-0',
            'E3SM-2-0',
            'E3SM-2-0-NARRM',
            'EC-Earth3',
            'EC-Earth3-AerChem',
            'EC-Earth3-CC',
            'EC-Earth3-Veg',
            'EC-Earth3-Veg-LR',
            'FGOALS-f3-L',
            'GFDL-CM4',
            'GFDL-ESM4',
            'HadGEM3-GC31-MM',
            'MPI-ESM1-2-HR',
            'MRI-ESM2-0',
            'NorESM2-MM',
            'SAM0-UNICON',
	    'TaiESM1',)

nc_name='pr'	    
unit_const=86400 # to convert to mm day^{-1}

# original netcdf names for latitude and longitude
tim_name='time'
lat_name='lat'
lon_name='lon'

python ITCZ_state_algorithm_absthresh_CMIP6_savenc.py 'TaiESM1,' $nc_name $tim_name $lat_name $lon_name $unit_const
#for item in "${my_array[@]}"; do
    #python ITCZ_state_algorithm_absthresh_CMIP6_savenc.py $item $nc_name $tim_name $lat_name $lon_name $unit_const

done	    
