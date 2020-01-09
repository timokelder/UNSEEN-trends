#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 14:36:43 2019

@author: timok
"""
#import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.stats import spearmanr
import pandas as pd
import xarray as xr
import os


#Domain for the West Coast (WC) and Svalbard (SV) domains
lats= [58,63]# SV[76,80]
lons=[4,7] #SV[8,30]

season=[8,9,10,11] #A+SON Autumn 
initialization_month=9

time_event=3 #number of day of the event

dirname=r'//home/timok/timok/SALIENSEAS/SEAS5'
#dirname=r'/home/timok/Documents/ensex/python'

ensemble_amount=25


season_start,season_end=str(2014)+'-'+str(season[0])+'-01',str(2014)+'-'+str(season[len(season)-1]+1)+'-01'  #Because the precipitation is cumulative, we take the 1st day of the month after the season as end date for decummulation

#Set directory
os.chdir(dirname)

#Open ensemble member 00
xr_0=xr.open_dataset(dirname+'/00/Arctic.SEAS5_sfc_00_20140801.nc')

#Select the region and reduce the domain season
xr_season=xr_0.sel(time=slice(season_start,season_end),lon=slice(lons[0],lons[1]),lat=slice(lats[0],lats[1])).drop(['SSTK','CI']) #drop the unnecessary variables
    
#open the netcdf with 200-year return quantiles to mask regions
xr_q200=xr.open_dataset(dirname+'/ensex/statistics/multiday/Quantile_ld2/Quantiles200.nc')

#The region is selected based on where the climatology (200-year values) are greater than a user-defined threshold
climatology_threshold=90 #35 #

#Select the regional averaged 3 day cumulative precipitation
xr_season_region=( xr_season['LSP'] #select the large scale precipitation variable
#.diff('time') #the LSP is cumulative precipitation -> convert this to daily precip
#.rolling(time=3).sum() #our target events are 3-day events -> convert daily to 3 day cumulative
.where(xr_q200['P_200']>climatology_threshold).mean(dim=['lat','lon']) #Select the regional average where the climatology is greater than 90 
*1000 #.max(dim='time') Select the maximum value for this season  
)

for i in range(1,ensemble_amount):
    xr_i=xr.open_dataset(dirname+'/'+"%02d" % i +'/Arctic.SEAS5_sfc_'+"%02d" % i +'_20140801.nc')
    xr_=xr_i.sel(time=slice(season_start,season_end),lon=slice(lons[0],lons[1]),lat=slice(lats[0],lats[1])).drop(['SSTK','CI']) #drop the unnecessary variables
    xr_season_region=xr.concat([xr_season_region,xr_['LSP'].where(xr_q200['P_200']>climatology_threshold).mean(dim=['lat','lon'])*1000],'ensemble') #.diff('time').rolling(time=time_event).sum().dropna('time')   .max(dim='time')


##Decumulate
#Open ensemble member 00
xr_1_0=xr.open_dataset(dirname+'/00/Arctic.SEAS5_sfc_00_20140801.nc')
xr_1_1=xr.open_dataset(dirname+'/01/Arctic.SEAS5_sfc_01_20140801.nc')

#Select the region and reduce the domain season
season_start,season_end=str(2014)+'-'+str(season[1])+'-01',str(2014)+'-'+str(season[len(season)-1]+1)+'-01'  #Because the precipitation is cumulative, we take the 1st day of the month after the season as end date for decummulation

xr_season_1_0=xr_1_0.sel(time=slice(season_start,season_end),lon=slice(lons[0],lons[1]),lat=slice(lats[0],lats[1])).drop(['SSTK','CI']) #drop the unnecessary variables
xr_season_1_1=xr_1_1.sel(time=slice(season_start,season_end),lon=slice(lons[0],lons[1]),lat=slice(lats[0],lats[1])).drop(['SSTK','CI']) #drop the unnecessary variables
    
#Select the regional averaged 3 day cumulative precipitation
xr_season_region_1_0=xr_season_1_0['LSP'].diff('time').rolling(time=3).sum().dropna('time').where(xr_q200['P_200']>climatology_threshold).mean(dim=['lat','lon'])*1000 #.max(dim='time') Select the maximum value for this season
xr_season_region_1_1=xr_season_1_1['LSP'].diff('time').rolling(time=3).sum().dropna('time').where(xr_q200['P_200']>climatology_threshold).mean(dim=['lat','lon'])*1000 #.max(dim='time') Select the maximum value for this season


season_region_series_1_0 = xr_season_region_1_0.to_series()
season_region_series_1_1 = xr_season_region_1_1.to_series()
    
    
###Set graphic settigns
font_graphs = {'family' : "Times New Roman",
    'weight' : 'normal',
    'size'   : 14}

plt.rc('font', **font_graphs)
    
days = mdates.DayLocator()   # every year
months = mdates.MonthLocator()  # every month
months_fmt = mdates.DateFormatter('%Y-%m')   
    
###Graphic

xr_season_region.name='SON-3DP ensemble forecasts (mm)'
f=plt.figure(figsize=(9, 5.2))
grid = plt.GridSpec(5, 1, wspace=0, hspace=0.1)
ax0 = f.add_subplot(grid[:-2,0])
xr_season_region.sel(ensemble=slice(2,24)).diff('time').rolling(time=3).sum().dropna('time').plot.line(x='time',linewidth=0.5,color='black',alpha=0.7,ax=ax0,add_legend=False) #aspect=2,size=3. Include (exclude) for decumulated (cumulative) precipitation: .diff('time').rolling(time=3).sum().dropna('time') 
xr_season_region.sel(ensemble=0).diff('time').rolling(time=3).sum().dropna('time').plot.line(x='time',linewidth=0.5,color='blue',ax=ax0,add_legend=False) #aspect=2,size=3
xr_season_region.sel(ensemble=1).diff('time').rolling(time=3).sum().dropna('time').plot.line(x='time',linewidth=0.8,color='orange',ax=ax0,add_legend=False) #aspect=2,size=3
ax0.yaxis.set_label_position("right")
ax0.yaxis.tick_right()
ax0.spines['top'].set_visible(False)
ax0.spines['left'].set_visible(False)
ax0.tick_params(axis='both',labelleft=False,labelbottom=False,labelright=True)
ax0.set_xlabel('')


ax1 = f.add_subplot(grid[-2,0], sharex=ax0)
plt.bar(season_region_series_1_0.index,season_region_series_1_0.values,color='blue')
ax1.tick_params(axis='both',labelleft=False,labelbottom=False,labelright=True)
ax1.set_xticks([], [])
ax1.yaxis.tick_right()
ax1.set_ylim(0,100)
ax1.spines['top'].set_visible(False)
ax1.spines['left'].set_visible(False)
#ax1.set_ylabel("Three-day")
#ax1.yaxis.set_label_position("right")

ax2 = f.add_subplot(grid[-1,0], sharex=ax0)
plt.bar(season_region_series_1_1.index,season_region_series_1_1.values,color='orange')
ax2.tick_params(axis='both', labelleft=False,labelright=True) #labelsize=18,labelbottom=False,
ax2.yaxis.tick_right()
ax2.spines['top'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.set_ylim(0,100)


ax2.xaxis.set_major_locator(months)
ax2.xaxis.set_major_formatter(months_fmt)
ax2.xaxis.set_minor_locator(days)
#ax.fmt_xdata = mdates.DateFormatter('%Y-%m')
#ax0.set_title('Cumulative')
#ax1.set_title('Three-day')

plt.xticks(rotation=0)
plt.margins(x=0)
plt.ylab='SON-3DP ensemble forecasts (mm)'
#f.text(0.97, 0.26, 'Three-day', ha='center', va='center', rotation='vertical')
f.text(1, 0.5, 'SON-3DP ensemble forecasts (mm)', ha='center', va='center', rotation='vertical')

#plt.tight_layout()
plt.savefig(dirname+'/ensex/statistics/multiday/plots/Ensemble_forecasts.svg',bbox_inches='tight')


## Artificial series over 35 years
Extremes=xr.open_dataset(dirname+'/ensex/Extremes/Extremes.nc')

MBR0=Extremes.sel(leadtime=2,ensemble=0).to_array().values.flatten()
MBR1=Extremes.sel(leadtime=2,ensemble=1).to_array().values.flatten()
MBR0_anomaly=(MBR0-MBR0.mean())/MBR0.std()
MBR1_anomaly=(MBR1-MBR1.mean())/MBR1.std()

f=plt.figure(figsize=(9, 3))#f=plt.figure(figsize=(20, 3))
grid = plt.GridSpec(2, 3, wspace=0.4, hspace=0.1)
ax0 = f.add_subplot(grid[0,:-1])

#ax0 = f.add_subplot(211)
Extremes.sel(leadtime=2,ensemble=0).to_array().plot(color='blue',marker='o',ax=ax0)
ax0.tick_params(axis='both',labelleft=False,labelbottom=False,labelright=True)# labelsize=18,
ax0.set_title('')
ax0.yaxis.tick_right()
ax1 = f.add_subplot(grid[1,:-1])
#ax1 = f.add_subplot(212)
Extremes.sel(leadtime=2,ensemble=1).to_array().plot(color='orange',marker='+',ax=ax1)
ax1.tick_params(axis='both', labelleft=False,labelbottom=True,labelright=True)
ax1.set_title('')
ax0.set_xlabel('')
ax1.set_xlabel('')
ax1.yaxis.tick_right()
ax2 = f.add_subplot(grid[:,-1])
plt.scatter(MBR0_anomaly,MBR1_anomaly,color='grey')
plt.plot([-2.5,2.5],[-2.5,2.5],color='black')
ax2.yaxis.tick_right()
ax2.tick_params(axis='both', labelleft=False,labelright=True) #labelsize=18,labelbottom=False,
ax2.tick_params(axis='x', colors='blue')
ax2.tick_params(axis='y', colors='orange')
ax2.text(-2.4,1.5,'r = '+"%1.2f" % spearmanr(MBR0_anomaly,MBR1_anomaly)[0])
#plt.show()

plt.savefig(dirname+'/ensex/statistics/multiday/plots/Correlate.svg',bbox_inches='tight')

