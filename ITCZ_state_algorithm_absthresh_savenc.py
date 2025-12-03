

# First load all necessary modules
import numpy as np 
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import sys

# folders where precipitation data lives
main_fld = ''

def main():
    if len(sys.argv) > 1:
        arguments = sys.argv[1:]
        print("Arguments received:", arguments)
        # Add your logic here to process the arguments
        for arg in arguments:
            print(f"Processing argument: {arg}")

        abs_thresh = float(arguments[0])
        acron = arguments[1]
        dataset_fld = arguments[2]
        fname    = arguments[3]
        nc_name  = arguments[4]
        tim_name = arguments[5]
        lat_name = arguments[6]
        lon_name = arguments[7]
        unit_const = float(arguments[8])
        region = arguments[9]
        region_str = arguments[10]
        lon0 = float(arguments[11])
        lon1 = float(arguments[12])
        
    else:
        print("No arguments provided.")
    return abs_thresh, acron, dataset_fld, fname, nc_name, tim_name, \
            lat_name, lon_name, unit_const, region, region_str, lon0, lon1

#print(yr0)
if __name__ == "__main__":
    [abs_thresh, acron, dataset_fld, fname, nc_name, tim_name, lat_name, \
            lon_name, unit_const, region, region_str, lon0, lon1] = main()

save_nc = True

#print(abs_thresh)
#print(acron)
#print(tim_name)
#print(unit_const)
#exit()

dri = main_fld+dataset_fld
fn = fname
data = xr.open_dataset(dri+fn)
print(acron)

# new names for latitude and longitude
tim_name_f = tim_name
lat_name_f = 'lat'
lon_name_f = 'lon'

#----threshold percentage of longitudes with dominant ITCZ state----#

# zonal running average
runavg_gpts = 5

# 30% implies 13/46 longitudes
# 35% implies 16/46 longitudes
# 40% implies 18/46 longitudes (~2000 km)
# 45% implies 20/46 longitudes
# 50% implies 23/46 longitudes (~2500 km)
# 58% implies 27/46 longitudes (~3000 km)
lon_thresh = 58

#----threshold percentile of regional precipitation (20S-20N by 45 deg wide domain)------#
# 85% implies 282/1886 grid points (~6 latitudes per longitude) are above threshold
# 90% implies 188/1886 grid points (~4 latitudes per longitude) are above threshold
# 95% implies 94/1886 grid points (~2 latitudes per longitude) are above threshold
prc_threshold = 90

#----------------------------------------------------------------#
##### THE USER MUST KNOW THE ORDER OF LATITUDE AND LONGITUDE #####
##### AND ORDER ALL VALUES BELOW CONSISTENT WITH THIS ORDER #####
#----------------------------------------------------------------#

# "tropical" latitudes (may need to be changed for different grids)
dlat = 1
lat0 = -20
lat1 = 20

# off-equatorial latitudes – must be positive only (may need to be changed for different grids)
lat0_offEQ = 2
lat1_offEQ = 20

# equatorial latitudes (may need to be changed for different grids)
lat0_EQ = -2+dlat
lat1_EQ = 2-dlat

#----------------------------------------------------------------#
#----------------------------------------------------------------#

# output directory and file
dro = dri
fn_output = acron+'_ITCZ_st_daily_'+region+'_abs_thresh_'+str(abs_thresh)+'mm'
print(dro+fn_output+'.nc')

# rename latitude and longitude and convert to mm/day
if (lat_name != lat_name_f) & (lon_name != lon_name_f):
    pri = data[nc_name].sortby(data[tim_name_f]).rename({lat_name: lat_name_f,lon_name: lon_name_f})*unit_const
else:
    pri = data[nc_name].sortby(data[tim_name_f])*unit_const

timei = pri[tim_name_f]
lati = pri[lat_name_f]
loni = pri[lon_name_f]
[nti,nlati,nloni] = pri.shape

# rename latitude and longitude and convert to mm/day
if (lat_name != lat_name_f) & (lon_name != lon_name_f):
    pri = data[nc_name].sortby(timei).rename({lat_name: lat_name_f,lon_name: lon_name_f})*unit_const
else:
    pri = data[nc_name].sortby(timei)*unit_const

# now enure longitudes are correct for region of interest
if (region == 'WestPac') & (loni.min() < 0):
   # longitudes are from -180 to 180; adjust pr data to be for longitudes from 0 to 360
   cond_lt0 = pri[lon_name_f] < 0
   cond_ge0 = pri[lon_name_f] >= 0

   lon_neg = pri[lon_name_f].where(cond_lt0, drop=True).values + 360
   lon_pos = pri[lon_name_f].where(cond_ge0, drop=True).values
   nlonn = len(lon_neg)
   nlonp = len(lon_pos)

   # final longitude array
   lonf = np.zeros((nloni))
   lonf[:nlonp] = lon_pos
   lonf[nlonp:] = lon_neg
   lon = xr.DataArray(lonf,coords=[lonf],dims=[lon_name_f])

   # adjust precip array
   pr_neg = pri.where(cond_lt0, drop=True).values
   pr_pos = pri.where(cond_ge0, drop=True).values

   pr_new = np.zeros((nti,nlati,nloni))
   pr_new[:,:,:nlonp] = pr_pos
   pr_new[:,:,nlonp:] = pr_neg
   pr_xr = xr.DataArray(pr_new,coords=[timei,lati,lonf],dims=[tim_name,lat_name_f,lon_name_f])
   pr = pr_xr.drop_duplicates(dim=lon_name_f)
elif (region == 'Atlantic') & (loni.max() > 180):
    # longitudes are from 0 to 360; adjust pr data to be for longitudes from -180 to 180
    pr_above180 = pr.where(pr[lon_name_f] > 180, drop=True)

    cond_lt180 = pri[lon_name_f] < 180
    cond_ge180 = pri[lon_name_f] >=180

    lon_lt180 = pri[lon_name_f].where(cond_lt180, drop=True).values
    lon_ge180 = pri[lon_name_f].where(cond_ge180, drop=True).values - 360
    nlonlt = len(lon_lt180)
    nlonge = len(lon_ge180)

    # final longitude array
    lonf = np.zeros((nloni))
    lonf[:nlonge] = lon_ge180
    lonf[nlonge:] = lon_lt180
    lon = xr.DataArray(lonf,coords=[lonf],dims=[lon_name_f])

    # adjust precip array
    pr_lt180 = pri.where(cond_lt180, drop=True).values
    pr_ge180 = pri.where(cond_ge180, drop=True).values

    pr_new = np.zeros((nti,nlati,nloni))
    pr_new[:,:,:nlonp] = pr_ge180
    pr_new[:,:,nlonp:] = pr_lt180
    pr = xr.DataArray(pr_new,coords=[timei,lati,lonf],dims=[tim_name,lat_name_f,lon_name_f])
else:
    pr = pri

time = timei
dataset = acron

# subset precip only over our region of interest
pr_region = pr.sel(lat=slice(lat0,lat1),lon=slice(lon0,lon1))

# global running average precipitation and subset running average precip only our region of interest
# (for now it is over 2 degrees longitude west and 2 degrees east of grid points, 5 points)
pr_region_run = ( pr.sel(lat=slice(lat0,lat1)).rolling(lon=runavg_gpts, center=True, min_periods=1).mean() ).sel(lon=slice(lon0,lon1))
pr_region_run



# subset precip data into NH, SH, and EQ regions

# NH precipitation
pr_NH = pr_region_run.sel(lat=slice(lat0_offEQ,lat1_offEQ))

# SH precipitation
pr_SH = pr_region_run.sel(lat=slice(-lat1_offEQ,-lat0_offEQ))

# EQ precipitation
if lat0_EQ == lat1_EQ: # not recommended to just use only one equatorial grid point
    pr_EQ = pr_region_run.sel(lat=lat0_EQ)
else:
    pr_EQ = pr_region_run.sel(lat=slice(lat0_EQ,lat1_EQ))

nlat = len(pr_region[lat_name_f])
nlon = len(pr_region[lon_name_f])
nt = len(pr_region)

# compute precipitation thresholds
pr_region_flat = pr_region_run.stack(z=(lat_name_f, lon_name_f)).reset_index('z')

pr_NH_flat = pr_NH.stack(z=(lat_name_f, lon_name_f)).reset_index('z')
pr_SH_flat = pr_SH.stack(z=(lat_name_f, lon_name_f)).reset_index('z')
pr_EQ_flat = pr_EQ.stack(z=(lat_name_f, lon_name_f)).reset_index('z')


# now take 90th percentile of precipitation
pr_tropics_low = abs_thresh#pr_region_run.quantile(prc_threshold/100., dim=(lat_name_f, lon_name_f)) #.9008
pr_low = pr_region_run.quantile(prc_threshold/100., dim=(lat_name_f, lon_name_f))
pr_NH_low      = pr_NH.quantile(prc_threshold/100., dim=(lat_name_f, lon_name_f))
pr_SH_low      = pr_SH.quantile(prc_threshold/100., dim=(lat_name_f, lon_name_f))

# normalized precipitation 
pr_NH_norm = pr_NH/pr_tropics_low
pr_SH_norm = pr_SH/pr_tropics_low
pr_EQ_norm = pr_EQ/pr_tropics_low

# drop all gridpoints below precip threshold
pr_NH_norm_thresh = pr_NH_norm.where(pr_NH_norm>=1).mean(dim=lat_name_f,skipna=True)
pr_SH_norm_thresh = pr_SH_norm.where(pr_SH_norm>=1).mean(dim=lat_name_f,skipna=True)
pr_EQ_norm_thresh = pr_EQ_norm.where(pr_EQ_norm>=1).mean(dim=lat_name_f,skipna=True)

# conditions for high precipitation
cond1 = pr_NH_norm_thresh >= 1
cond2 = pr_SH_norm_thresh >= 1
cond3 = pr_EQ_norm_thresh >= 1

# conditions for low precipitation
cond_a1 = np.isnan(pr_NH_norm_thresh) == True
cond_a2 = np.isnan(pr_SH_norm_thresh) == True
cond_a3 = np.isnan(pr_EQ_norm_thresh) == True

#----------Northern ITCZ------------

# Condition 1: large NH, small EQ, and small SH
tf_nITCZ = np.where( (cond1 & cond_a2 & cond_a3),1,0 )
count_nITCZ = np.sum(tf_nITCZ,axis=1)
prc_nITCZ = count_nITCZ/nlon*100.

#----------Northern ITCZ------------


#----------Southern ITCZ------------

# Condition 2: large SH, small EQ, and small NH
tf_sITCZ = np.where( (cond2 & cond_a1 & cond_a3),1,0 )
count_sITCZ = np.sum(tf_sITCZ,axis=1)
prc_sITCZ = count_sITCZ/nlon*100.

#----------Southern ITCZ------------


#----------Double ITCZ--------------

# Condition 3: large SH, small EQ, and large NH
tf_dITCZ = np.where( (cond1 & cond2 & cond_a3),1,0 )
count_dITCZ = np.sum(tf_dITCZ,axis=1)
prc_dITCZ = count_dITCZ/nlon*100.

#----------Double ITCZ--------------


#----------Thin Equatorial ITCZ----------

# Condition 4: small SH, large EQ, and small NH
tf_eITCZ = np.where( (cond3 & cond_a1 & cond_a2),1,0 )
count_eITCZ = np.sum(tf_eITCZ,axis=1)
prc_eITCZ = count_eITCZ/nlon*100.

#----------Thin Equatorial ITCZ----------

#----------Wide Equatorial ITCZ (Gets combined with eITCZ)----------

# Condition 5: large NH, large EQ, and small SH 
tf_wITCZ_NH = np.where( (cond1 & cond_a2 & cond3),1,0 )
count_wITCZ_NH = np.sum(tf_wITCZ_NH,axis=1)
prc_wITCZ_NH = count_wITCZ_NH/nlon*100.

# Condition 6: large SH, large EQ, and small NH 
tf_wITCZ_SH = np.where( (cond2 & cond_a1 & cond3),1,0 )
count_wITCZ_SH = np.sum(tf_wITCZ_SH,axis=1)
prc_wITCZ_SH = count_wITCZ_SH/nlon*100.

# Condition 7: large SH, large EQ, and large NH 
tf_wITCZ_EQ = np.where( (cond1 & cond2 & cond3),1,0 )
count_wITCZ_EQ = np.sum(tf_wITCZ_EQ,axis=1)
prc_wITCZ_EQ = count_wITCZ_EQ/nlon*100.

count_wITCZ = count_wITCZ_NH + count_wITCZ_SH + count_wITCZ_EQ
prc_wITCZ = count_wITCZ/nlon*100.

#----------Wide Equatorial ITCZ----------

#----------Absent ITCZ----------

# Condition 8: small SH, small EQ, and small NH
tf_aITCZ = np.where( (cond_a1 & cond_a2 & cond_a3),1,0 )
count_aITCZ = np.sum(tf_aITCZ,axis=1)
prc_aITCZ = count_aITCZ/nlon*100.

#----------Absent ITCZ----------

# save a longitude based ITCZ state array for plotting purposes
tf_dITCZ_big = tf_dITCZ*1
tf_nITCZ_big = tf_nITCZ*2
tf_sITCZ_big = tf_sITCZ*3
tf_aITCZ_big = tf_aITCZ*4
tf_eITCZ_big = (tf_eITCZ + tf_wITCZ_NH + tf_wITCZ_SH + tf_wITCZ_EQ)*5

tf_ITCZ_big = tf_dITCZ_big + tf_nITCZ_big + tf_sITCZ_big + tf_aITCZ_big + tf_eITCZ_big
tf_ITCZ_big_xr = xr.DataArray(tf_ITCZ_big, coords=[time,pr_region[lon_name_f]], dims=[tim_name,lon_name_f])

tf_ITCZ_big_xr.name = 'itcz_state_lon'
tf_ITCZ_big_xr.attrs['long_name'] = dataset+' daily ITCZ state as a function of longitude'
tf_ITCZ_big_xr.attrs['ITCZ state numbers'] = '1=dITCZ, 2=nITCZ, 3=sITCZ, 4=aITCZ, 5=eITCZ'

prc_nITCZ_xr = xr.DataArray(prc_nITCZ, coords=[time], dims=[tim_name])
prc_sITCZ_xr = xr.DataArray(prc_sITCZ, coords=[time], dims=[tim_name])
prc_dITCZ_xr = xr.DataArray(prc_dITCZ, coords=[time], dims=[tim_name])
prc_eITCZ_xr = xr.DataArray(prc_eITCZ, coords=[time], dims=[tim_name])
prc_wITCZ_xr = xr.DataArray(prc_wITCZ, coords=[time], dims=[tim_name])
prc_aITCZ_xr = xr.DataArray(prc_aITCZ, coords=[time], dims=[tim_name])

itcz_state_num = np.arange(1,7,1)
itcz_state_num = xr.DataArray(itcz_state_num, coords=[itcz_state_num], dims=['itcz_state_number'])

prc_all = xr.DataArray(np.zeros((6,nt)), coords=[itcz_state_num,time], dims=['itcz_state_number',tim_name])
prc_all[0,:] = prc_dITCZ
prc_all[1,:] = prc_nITCZ
prc_all[2,:] = prc_sITCZ
prc_all[3,:] = prc_aITCZ
prc_all[4,:] = prc_eITCZ
prc_all[5,:] = prc_wITCZ

#--------------THIS PART IS SOMEWHAT SLOW-----------------#
# Choose dominant ITCZ state
itcz_state_new = np.zeros((nt))
for tt in range(0,nt,1):
    prc_max = prc_all[:,tt].max()
    #print(prc_max.values)
    itcz_state_max = itcz_state_num.where(prc_all[:,tt]==prc_max.values,drop=True)
    itcz_state_new[tt] = itcz_state_max[0].values
#--------------THIS PART IS SOMEWHAT SLOW-----------------#

# make itcz state array an xarray
itcz_state_new_xr = xr.DataArray(itcz_state_new, coords=[time], dims=[tim_name])

prc_max_xr = prc_all.max(dim=('itcz_state_number')) # over ITCZ state
#prc_max_xr

# adding wITCZ_NH and prc_wITCZ_EQ to dITCZ
prc_all_ndI = prc_dITCZ + prc_wITCZ_NH + prc_wITCZ_EQ 

# adding wITCZ_SH and prc_wITCZ_EQ to dITCZ
prc_all_sdI = prc_dITCZ + prc_wITCZ_SH + prc_wITCZ_EQ 

# If an ITCZ is nITCZ but also has some longitudes with dITCZ or wITCZ_NH,
# add those longitudes to nITCZ percentage
prc_max_final = np.where(itcz_state_new_xr==2, prc_max_xr + prc_all_ndI, prc_max_xr) # add dITCZ and wITCZ_NH to nITCZ

# If an ITCZ is sITCZ but also has some longitudes with dITCZ or wITCZ_SH,
# add those longitudes to sITCZ percentage
prc_max_final = np.where(itcz_state_new_xr==3, prc_max_xr + prc_all_sdI, prc_max_final) # add dITCZ and wITCZ_SH to sITCZ

# Now convert updated percentages to an xarray
prc_max_final_xr = xr.DataArray(prc_max_final, coords=[time], dims=[tim_name])

# Any days below longitude percent threshold will be assigned aITCZ
itcz_state_new_xr = xr.where(prc_max_final_xr<lon_thresh,4,itcz_state_new_xr)

itcz_state_new_xr
itcz_state_new_xr.name = 'itcz_state'
itcz_state_new_xr.attrs['long_name'] = dataset+' daily ITCZ state'
itcz_state_new_xr.attrs['region'] = region_str
itcz_state_new_xr.attrs['method'] = 'Algorithm developed for daily averaged precipitation rate'


itcz_state_new_xr_final = itcz_state_new_xr*0 + np.where(itcz_state_new_xr.values==6,5,itcz_state_new_xr.values)
itcz_state_new_xr_final

itcz_state_new_xr_final.name = 'itcz_state'
itcz_state_new_xr_final.attrs['long_name'] = dataset+' daily ITCZ state'
itcz_state_new_xr_final.attrs['method'] = 'Algorithm developed for daily averaged precipitation rate'
itcz_state_new_xr_final.attrs['region'] = region_str
itcz_state_new_xr_final.attrs['precip_threshold'] = str(abs_thresh)+' mm/day'
itcz_state_new_xr_final.attrs['longitude_threshold'] = str(lon_thresh)+' percent of gridpoints must agree on ITCZ state'
itcz_state_new_xr_final.attrs['ITCZ state numbers'] = '1=dITCZ, 2=nITCZ, 3=sITCZ, 4=aITCZ, 5=eITCZ'

pr_tropics_low_final = pr_low
pr_tropics_low_final.name = 'pr_threshold'
pr_tropics_low_final.attrs['long_name'] = dataset+' precipitation threshold'
pr_tropics_low_final.attrs['region'] = region_str
pr_tropics_low_final.attrs['precip_threshold'] = str(prc_threshold)+'th percentile of regional precip'

pr_SH_low_final = pr_SH_low
pr_SH_low_final.name = 'pr_threshold_SH'
pr_SH_low_final.attrs['long_name'] = dataset+' precipitation threshold in SH (2-20S)'
pr_SH_low_final.attrs['region'] = region_str
pr_SH_low_final.attrs['precip_threshold'] = str(prc_threshold)+'th percentile of regional precip'

pr_NH_low_final = pr_NH_low
pr_NH_low_final.name = 'pr_threshold_NH'
pr_NH_low_final.attrs['long_name'] = dataset+' precipitation threshold in NH (2-20N)'
pr_NH_low_final.attrs['region'] = region_str
pr_NH_low_final.attrs['precip_threshold'] = str(prc_threshold)+'th percentile of regional precip'

dataset_final = xr.merge([itcz_state_new_xr_final,pr_tropics_low_final,pr_SH_low_final,pr_NH_low_final,tf_ITCZ_big_xr])
print(dataset_final)

# now composite average precip and compute % of days for specific ITCZ states
# set up conditional statements 
cond_dITCZ = (itcz_state_new_xr_final == 1)
cond_nITCZ = (itcz_state_new_xr_final == 2)
cond_sITCZ = (itcz_state_new_xr_final == 3)
cond_aITCZ = (itcz_state_new_xr_final == 4)
cond_eITCZ = (itcz_state_new_xr_final == 5)

# compute percentages for each ITCZ state
prc_dITCZ = round(len(itcz_state_new_xr_final.where((cond_dITCZ), drop=True))/(len(time))*100,1)
prc_nITCZ = round(len(itcz_state_new_xr_final.where((cond_nITCZ), drop=True))/(len(time))*100,1)
prc_sITCZ = round(len(itcz_state_new_xr_final.where((cond_sITCZ), drop=True))/(len(time))*100,1)
prc_aITCZ = round(len(itcz_state_new_xr_final.where((cond_aITCZ), drop=True))/(len(time))*100,1)
prc_eITCZ = round(len(itcz_state_new_xr_final.where((cond_eITCZ), drop=True))/(len(time))*100,1)

print('Double     ITCZ happens '+str(prc_dITCZ)+'% of the time.')
print('Northern   ITCZ happens '+str(prc_nITCZ)+'% of the time.')
print('Southern   ITCZ happens '+str(prc_sITCZ)+'% of the time.')
print('Absent     ITCZ happens '+str(prc_aITCZ)+'% of the time.')
print('Equatorial ITCZ happens '+str(prc_eITCZ)+'% of the time.')

# save netcdf ITCZ states
if save_nc == True:
    dataset_final.to_netcdf(path=dro+fn_output+'.nc')
