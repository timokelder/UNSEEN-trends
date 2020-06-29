#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 15:38:59 2019

@author: timok
"""
import os
dirname=r'/Users/Timo/OneDrive - Loughborough University/GitHub/UNSEEN-trends'
os.chdir(dirname)

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import gridspec
from mpl_toolkits.basemap import Basemap

# import netCDF4 as nc
import pandas as pd


plt.rcParams["font.family"] = "cursive" 
plt.rcParams['font.size'] = 7

 
m = Basemap(llcrnrlon=2.,urcrnrlon=20.,llcrnrlat=57.2,urcrnrlat=70.,
             resolution='i', projection='lcc', lat_1=57., lon_0=5.)
    
m_sv = Basemap(llcrnrlon=7.,urcrnrlon=35.,llcrnrlat=75.,urcrnrlat=81.,
             resolution='i', projection='lcc', lat_1=65., lon_0=5.)

cm_data_oslo = np.loadtxt("ScientificColourMaps4/oslo/oslo.txt") #vik#imola or oslo
Oslo=LinearSegmentedColormap.from_list('CBname', cm_data_oslo)[::-1]
statistics=pd.read_table('/home/timok/ensex/statistics/multiday/Quantile_ld2/statistics_quantile200.txt',sep=',')#[,2:4]

fig = plt.figure(linewidth=5.5)
ax = fig.add_axes([0.044,0.05,0.9,0.9])
statistics_threshold =np.array(statistics['Quantile 200-yr'])
im1 = m.pcolormesh(np.array(statistics['lon']).reshape(100,141),np.array(statistics['lat']).reshape(100,141),np.array(statistics['Quantile 200-yr']).reshape(100,141),shading='flat',cmap=Oslo,latlon=True)#vmin=vmin_vmax[k][0],vmax=vmin_vmax[k][1],
m.contour(np.array(statistics['lon']).reshape(100,141),np.array(statistics['lat']).reshape(100,141),np.array(statistics['Quantile 200-yr']).reshape(100,141),levels=[90],colors='white',latlon=True)#vmin=vmin_vmax[k][0],vmax=vmin_vmax[k][1],
x1,y1=m([4,4,7,7],[58,63,63,58])
#m.scatter(x1,y1,s=1, c='k', marker='.', alpha=.5) 
xy = zip(x1,y1)
poly = matplotlib.patches.Polygon( xy, edgecolor='red', facecolor='None')
plt.gca().add_patch(poly)
#ax.set_title('SEAS5 RV2 (mm/3days)')
cb = m.colorbar(im1,"right", size="5%", pad="2%")#m.readshapefile('gadm36_NOR_shp/gadm36_NOR_0',name='gadm36_NOR_0',drawbounds=True,color='black')
m.drawcoastlines()
m.drawcountries()
#m.pcolor(lons_merge_grid[mask_WC_merge].reshape(11,9),lats_merge_grid[mask_WC_merge].reshape(11,9),RV20_seas5[mask_WC_merge].reshape(11,9),latlon=True, hatch='///', alpha=0.5, vmin=80, vmax=3000)
#m.contour(lons_merge_grid[mask_WC_merge].reshape(11,9),lats_merge_grid[mask_WC_merge].reshape(11,9),RV20_seas5[mask_WC_merge].reshape(11,9),levels=[0,80,300],latlon=True)
#m.contourf(lons_merge_grid[mask_WC_merge].reshape(11,9),lats_merge_grid[mask_WC_merge].reshape(11,9),RV20_seas5[mask_WC_merge].reshape(11,9),levels=[0,50,280],edgecolor='blue',facecolor='blue',color='blue',latlon=True)

parallels = np.arange(0.,90,5.)
par=m.drawparallels(parallels,labels=[1,0,0,0],fontsize=14)#,xoffset= is not working?
x0,x1 = ax.get_xlim()
w = x1-x0
for key, (lines,texts) in par.items():
    for text in texts:
        x,y = text.get_position()
        text.set_position((x0-0.1*w,y))

meridians = np.arange(-180.,180.,10.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=14)                     

#ax.set_title('RV 20')
ax2=plt.axes([0.26,0.54,0.3,0.4])
im2 = m_sv.pcolormesh(np.array(statistics['lon']).reshape(100,141),np.array(statistics['lat']).reshape(100,141),np.array(statistics['Quantile 200-yr']).reshape(100,141),shading='flat',vmin=10,vmax=60,cmap=Oslo,latlon=True)
m_sv.contour(np.array(statistics['lon']).reshape(100,141),np.array(statistics['lat']).reshape(100,141),np.array(statistics['Quantile 200-yr']).reshape(100,141),colors='white',levels=[35],latlon=True)#vmin=vmin_vmax[k][0],vmax=vmin_vmax[k][1],
x1_s,y1_s=m_sv([8,8,30,30],[76,80.25,80.25,76])
#m.scatter(x1,y1,s=1, c='k', marker='.', alpha=.5) 
xy_s = zip(x1_s,y1_s)
poly_s = matplotlib.patches.Polygon( xy_s, edgecolor='red', facecolor='None')
plt.gca().add_patch(poly_s)
cb = m_sv.colorbar(im2,"left", size="5%", pad="18%")#m.readshapefile('gadm36_NOR_shp/gadm36_NOR_0',name='gadm36_NOR_0',drawbounds=True,color='black')
#cb.ax.yaxis.set_tick_params(color='white')
plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color='white')
m_sv.drawcoastlines()
m_sv.drawcountries()
# draw parallels and meridians.
parallels = np.arange(0.,90,5.)
m_sv.drawparallels(parallels,labels=[0,0,0,0])#[0,1,0,0]
meridians = np.arange(-180.,180.,10.)
m_sv.drawmeridians(meridians,labels=[0,0,0,0])#[0,0,0,1]

size=fig.get_size_inches()
fig.set_size_inches(size*2)

fig.savefig('/home/timok/Documents/ensex/python/graphs/3day/Domain_Quantile200.png',bbox_inches='tight')

