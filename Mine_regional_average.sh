 /bin/bash
#$ -l h_rt=2:00:00
#$ -q ded-parallelx.q
#$ -l h_vmem=3G
#$ -t 1-35
#$ -o //home/timok/timok/ensex/out/out_$JOB_NAME.$JOB_ID.$HOSTNAME.$TASK_ID
#$ -e //home/timok/timok/ensex/err/err_$JOB_NAME.$JOB_ID.$HOSTNAME.$TASK_ID
#$ -cwd

echo "Got $NSLOTS slots for job $SGE_TASK_ID."

cat > "//home/timok/timok/ensex/merge/mine_SEAS5_region""$SGE_TASK_ID"".py" << EOF
######################################################################
######################################################################

import xarray as xr
import numpy as np
import os
import glob

#######Start######

#Domain for the West Coast (WC) and Svalbard (SV) domains
lats= [76,80] #58,63]# SV[76,80]
lons=[8,30] #4,7] #SV[8,30]

#Define the season and event duration
#Watch out with the season selection. In the script I use month-6, so month 1-6 should be 13-18. 
#Not sure if 13-18 work correctly, needs to be checked when using this.
season=[9,10,11] #SON Autumn 

time_event=3 #number of day of the event

dirname=r'//home/timok/timok/SALIENSEAS/SEAS5'

#def main(input_args):  ##Make more neat and possible for anyone to run
#    alles hier
#    
#if __name__ == '__main__':
#    parser = argparse.ArgumentParser(description='Extract the regional averaged seasonal extremes from the hindcast')
#    parser.add_argument('--lats', required=True,
#                        help='The latitudes over which the ')
#    args = parser.parse_args()
#    main(args)
    
#amount of ensembles
ensemble_amount=25

#Amount of initialization dates (in lead times)
lds=np.arange(2,2+7-len(season),1) #The files run for 7 monts: the longest lead times still forecasting the end of the target season depends on the length of the season selected. First lead time is removed 
initialization_months=season[len(season)-1]-lds-1 #

#Create a string for slicing the season
season_start,season_end=str($SGE_TASK_ID+1980)+'-'+str(season[0])+'-01',str($SGE_TASK_ID+1980)+'-'+str(season[len(season)-1]+1)+'-01'  #Because the precipitation is cumulative, we take the 1st day of the month after the season as end date for decummulation

#Set directory
os.chdir(dirname)


#open the 00 member forecast with Xarray for the different initialization months (in lead times)

######## NOT NEAT  ##########
xr_2=xr.open_dataset(dirname+'/00/Arctic.SEAS5_sfc_00_'+str($SGE_TASK_ID+1980)+"%02d" % initialization_months[0]+'01.nc')
xr_3=xr.open_dataset(dirname+'/00/Arctic.SEAS5_sfc_00_'+str($SGE_TASK_ID+1980)+"%02d" % initialization_months[1]+'01.nc')
xr_4=xr.open_dataset(dirname+'/00/Arctic.SEAS5_sfc_00_'+str($SGE_TASK_ID+1980)+"%02d" % initialization_months[2]+'01.nc')
xr_5=xr.open_dataset(dirname+'/00/Arctic.SEAS5_sfc_00_'+str($SGE_TASK_ID+1980)+"%02d" % initialization_months[3]+'01.nc')


#Select the domain and the season for the forecast
xr2=( xr_2.sel(time=slice(season_start,season_end),  #Select season
             lon=slice(lons[0],lons[1]),lat=slice(lats[0],lats[1])) #Select domain
.drop(['SSTK','CI']) #drop the unnecessary variables
)

#Same for the other lead times
xr3=xr_3.sel(time=slice(season_start,season_end),lon=slice(lons[0],lons[1]),lat=slice(lats[0],lats[1])).drop(['SSTK','CI']) #drop the unnecessary variables
xr4=xr_4.sel(time=slice(season_start,season_end),lon=slice(lons[0],lons[1]),lat=slice(lats[0],lats[1])).drop(['SSTK','CI']) #drop the unnecessary variables
xr5=xr_5.sel(time=slice(season_start,season_end),lon=slice(lons[0],lons[1]),lat=slice(lats[0],lats[1])).drop(['SSTK','CI']) #drop the unnecessary variables

    
#open the netcdf with 200-year return quantiles to mask regions
xr_q200=xr.open_dataset(dirname+'/ensex/statistics/multiday/Quantile_ld2/Quantiles200.nc')

#The region is selected based on where the climatology (200-year values) are greater than a user-defined threshold
climatology_threshold=35 #90

#Select the regional averaged 3 day cumulative precipitation
xr_ens2=( xr2['LSP'] #select the large scale precipitation variable
.diff('time') #the LSP is cumulative precipitation -> convert this to daily precip
.rolling(time=3).sum().dropna('time') #our target events are 3-day events -> convert daily to 3 day cumulative. Remove the NA that result from the rolling window
.where(xr_q200['P_200']>climatology_threshold).mean(dim=['lat','lon']) #Select the regional average where the climatology is greater than 90 
.max(dim='time')*1000 #Select the maximum value for this season
)

#Same for the other lead times
xr_ens3=xr3['LSP'].diff('time').rolling(time=3).sum().dropna('time').where(xr_q200['P_200']>climatology_threshold).mean(dim=['lat','lon']).max(dim='time')*1000
xr_ens4=xr4['LSP'].diff('time').rolling(time=3).sum().dropna('time').where(xr_q200['P_200']>climatology_threshold).mean(dim=['lat','lon']).max(dim='time')*1000
xr_ens5=xr5['LSP'].diff('time').rolling(time=3).sum().dropna('time').where(xr_q200['P_200']>climatology_threshold).mean(dim=['lat','lon']).max(dim='time')*1000


for i in range(1,ensemble_amount):
#    print(i)
    xr2_i=xr.open_dataset(dirname+'/'+"%02d" % i +'/Arctic.SEAS5_sfc_'+"%02d" % i +'_'+str($SGE_TASK_ID+1980)+"%02d" % initialization_months[0]+'01.nc')
    xr3_i=xr.open_dataset(dirname+'/'+"%02d" % i +'/Arctic.SEAS5_sfc_'+"%02d" % i +'_'+str($SGE_TASK_ID+1980)+"%02d" % initialization_months[1]+'01.nc')
    xr4_i=xr.open_dataset(dirname+'/'+"%02d" % i +'/Arctic.SEAS5_sfc_'+"%02d" % i +'_'+str($SGE_TASK_ID+1980)+"%02d" % initialization_months[2]+'01.nc')
    xr5_i=xr.open_dataset(dirname+'/'+"%02d" % i +'/Arctic.SEAS5_sfc_'+"%02d" % i +'_'+str($SGE_TASK_ID+1980)+"%02d" % initialization_months[3]+'01.nc')
    xr2=xr2_i.sel(time=slice(season_start,season_end),lon=slice(lons[0],lons[1]),lat=slice(lats[0],lats[1])).drop(['SSTK','CI']) #drop the unnecessary variables
    xr3=xr3_i.sel(time=slice(season_start,season_end),lon=slice(lons[0],lons[1]),lat=slice(lats[0],lats[1])).drop(['SSTK','CI']) #drop the unnecessary variables
    xr4=xr4_i.sel(time=slice(season_start,season_end),lon=slice(lons[0],lons[1]),lat=slice(lats[0],lats[1])).drop(['SSTK','CI']) #drop the unnecessary variables
    xr5=xr5_i.sel(time=slice(season_start,season_end),lon=slice(lons[0],lons[1]),lat=slice(lats[0],lats[1])).drop(['SSTK','CI']) #drop the unnecessary variables
    xr_ens2=xr.concat([xr_ens2,xr2['LSP'].diff('time').rolling(time=time_event).sum().dropna('time').where(xr_q200['P_200']>climatology_threshold).mean(dim=['lat','lon']).max(dim='time')*1000],'ensemble')
    xr_ens3=xr.concat([xr_ens3,xr3['LSP'].diff('time').rolling(time=time_event).sum().dropna('time').where(xr_q200['P_200']>climatology_threshold).mean(dim=['lat','lon']).max(dim='time')*1000],'ensemble')
    xr_ens4=xr.concat([xr_ens4,xr4['LSP'].diff('time').rolling(time=time_event).sum().dropna('time').where(xr_q200['P_200']>climatology_threshold).mean(dim=['lat','lon']).max(dim='time')*1000],'ensemble')
    xr_ens5=xr.concat([xr_ens5,xr5['LSP'].diff('time').rolling(time=time_event).sum().dropna('time').where(xr_q200['P_200']>climatology_threshold).mean(dim=['lat','lon']).max(dim='time')*1000],'ensemble')

 
xr_ens_lds=xr.concat([xr_ens2,xr_ens3,xr_ens4,xr_ens5],'leadtime')
xr_ens_lds=xr_ens_lds.assign_coords(leadtime=lds).assign_coords(ensemble=range(ensemble_amount))

xr_ens_lds.to_netcdf('//home/timok/timok/SALIENSEAS/SEAS5/ensex/Extremes/Ens_lds_SV'+str($SGE_TASK_ID)+'.nc', mode='w')


EOF

python "//home/timok/timok/ensex/merge/mine_SEAS5_region""$SGE_TASK_ID"".py"

