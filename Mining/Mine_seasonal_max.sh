 /bin/bash
#$ -l h_rt=2:00:00
#$ -q ded-parallelx.q
#$ -l h_vmem=3G
#$ -t 1-35
#$ -o //home/timok/timok/ensex/out/out_$JOB_NAME.$JOB_ID.$HOSTNAME.$TASK_ID
#$ -e //home/timok/timok/ensex/err/err_$JOB_NAME.$JOB_ID.$HOSTNAME.$TASK_ID
#$ -cwd

echo "Got $NSLOTS slots for job $SGE_TASK_ID."

cat > "//home/timok/timok/ensex/merge/mine_SEAS5""$SGE_TASK_ID"".py" << EOF
######################################################################
######################################################################
import netCDF4 as nc
import os
import numpy as np
import datetime
import calendar
import glob
from dateutil.relativedelta import relativedelta

#######Start######
lat_start= 57#latitude and longitude from https://www.gps-coordinates.net/
lat_end=82 #latitudes.units = 'degree_north'
lon_start=0 #longitudes.units = 'degree_east' 
lon_end=35
ensemble_amount=25 #will be 25

#Watch out with the season selection. In the script I use month-6, so month 1-6 should be 13-18. 
#Not sure if 13-18 work correctly, needs to be checked when using this.
season=[9,10,11] #SON Autumn 
n_event=3 #number of day of the event

dirname=r'//home/timok/timok/SALIENSEAS/SEAS5'
os.chdir(dirname)

#extract all the files that forecast the target season. In this case, select 35 years, 25 members, 4 lead times: 3500 files. We split the analysis over the years, so just select 25 members and 4 lead times. 
#The files run for 7 monts: end of the target season -6 months is the first initialization date that forecasts the target season (-6 instead of 7 because it is initialized the first of the month and we want it until the end of the month of the target season)
#The   
files = glob.glob('//home/timok/timok/SALIENSEAS/SEAS5/00/Arctic.SEAS5_sfc_00_'+str($SGE_TASK_ID+1980)+"%02d" % (season[len(season)-1]-6)+'01.nc') #read in the first file for this year. In this case 11(November)-6 = may is the first initialization month
for i in (range(7-1-len(season))):  # we further want the initialization dates up till the first month of the season. There are 7 initialization dates minus the amount of months in the target season. In this case 4 initialization dates, we need 3 more.  
    print ("%02d" % (season[len(season)-1]-6+i+1)+'01.nc')
    files += glob.glob('//home/timok/timok/SALIENSEAS/SEAS5/00/Arctic.SEAS5_sfc_00_'+str($SGE_TASK_ID+1980)+"%02d" % (season[len(season)-1]-6+i+1)+'01.nc')
        

for j in range(1,ensemble_amount):
    for i in (range(7-len(season))):  # Now select the files for all ensembles. In this case 4 initialization dates, 25 members
        print ("%02d" % (season[len(season)-1]-6+i)+'01.nc')
        files += glob.glob('//home/timok/timok/SALIENSEAS/SEAS5/'+"%02d" % j +'/Arctic.SEAS5_sfc_'+"%02d" % j +'_'+str($SGE_TASK_ID+1980)+"%02d" % (season[len(season)-1]-6+i)+'01.nc')
        


infile=nc.Dataset(files[0],'r')
time = infile.variables['time']
date_start = datetime.datetime ($SGE_TASK_ID+1980,season[0],1)
date_end = datetime.datetime ($SGE_TASK_ID+1980,1+season[len(season)-1],1) #1+last month to include the entire month. Build in If season contains months>12: season[len(season)-1] - 12 (to get the right month), 1980+1 (to select the next year)
index_start = nc.date2index(date_start, time)
index_end = nc.date2index(date_end, time)

T2M = infile.variables['T2M'][index_start:index_end,:,:] # [time x lat x lon]
# Differences with a step of 1 along the 'time' axis (0) 
ntim, nlat, nlon = np.shape(T2M)
lat_in = infile.variables['lat'][:]
lon_in = infile.variables['lon'][:]
nlat_domain = len(infile.variables['lat'][(lat_in >= lat_start) & (lat_in <= lat_end)])
nlon_domain = len(infile.variables['lon'][(lon_in >= lon_start) & (lon_in <= lon_end)])
time_end = 1 #ntim if you want to store the entire season
time_start=0

##### Just for testing
#I do not use pandas for increased speed
#LSP = infile.variables['LSP'][index_start:index_end+1,(lat_in >= lat_start) & (lat_in <= lat_end),(lon_in >= lon_start) & (lon_in <= lon_end)]
#LSP_diff = LSP[n_event:,:,:] - LSP[:-n_event,:,:] # LSP is cumulative precipitation. Take the 3-day precip
#LSP_diffmax=np.max(LSP_diff, axis=0)*1000 #Take maximum 3day event in the season, convert m to mm
##np.diff(LSP, n=1, axis=0) for check, doesnt use window


ncfile_out = nc.Dataset('//home/timok/timok/SALIENSEAS/SEAS5/ensex/Aggregated_months/seasonSON'+str(date_start.year)+'.nc', 'w')
ncfile_out.createDimension('time', None) #per ensemble, there are 6 different initialization dates
ncfile_out.createDimension('lat', nlat_domain)
ncfile_out.createDimension('lon', nlon_domain) 

LSP_out = ncfile_out.createVariable('LSP', 'f4', ('time', 'lat', 'lon',))
#T2M_out = ncfile_out.createVariable('T2M', 'f4', ('time', 'lat', 'lon',))
latitudes = ncfile_out.createVariable('lat', np.float32, ('lat',))
longitudes = ncfile_out.createVariable('lon', np.float32, ('lon',))
latitudes[:] = infile.variables['lat'][(lat_in >= lat_start) & (lat_in <= lat_end)]
longitudes[:] = infile.variables['lon'][(lon_in >= lon_start) & (lon_in <= lon_end)]
infile.close()   

latitudes.units = 'degree_north'
longitudes.units = 'degree_east' 

for file in files:
#    print(file)
    ncfile=nc.Dataset(file,'r') 
    time = ncfile.variables['time']
    index_start = nc.date2index(date_start, time)
    index_end = nc.date2index(date_end, time)
    LSP = ncfile.variables['LSP'][index_start:index_end+1,(lat_in >= lat_start) & (lat_in <= lat_end),(lon_in >= lon_start) & (lon_in <= lon_end)] # [time x lat x lon] for decummulation, one extra is needed
    LSP_diff = LSP[n_event:,:,:] - LSP[:-n_event,:,:] # LSP is cumulative precipitation. Take the 3-day precip
    LSP_diffmax=np.max(LSP_diff, axis=0)*1000 #Take maximum 3day event in the season, convert m to mm
#    T2M = ncfile.variables['T2M'][index_start:index_end,(lat_in >= lat_start) & (lat_in <= lat_end),(lon_in >= lon_start) & (lon_in <= lon_end)] # [time x lat x lon]
    ncfile.close()
    # Write out the new variable to a new file     
    #print file
    LSP_out[time_start:time_end,:,:] = LSP_diffmax
#    T2M_out[time_start:time_end,:,:] = T2M[:,:,:]
    time_start+=1 #ntim if you want to store the entire season
    time_end+=1

    
ncfile_out.close()


EOF

python "//home/timok/timok/ensex/merge/mine_SEAS5""$SGE_TASK_ID"".py"

