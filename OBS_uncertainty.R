import xarray as xr
import numpy as np
import os
import glob

#######Start######

#Domain for the West Coast (WC) and Svalbard (SV) domains
lats= [58,63]# SV[76,80]
lons=[4,7] #SV[8,30]

#Define the season and event duration
#Watch out with the season selection. In the script I use month-6, so month 1-6 should be 13-18. 
#Not sure if 13-18 work correctly, needs to be checked when using this.
season=[9,10,11] #SON Autumn 

time_event=3 #number of day of the event

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