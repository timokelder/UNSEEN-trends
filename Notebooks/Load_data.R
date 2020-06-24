library('ncdf4')
library('ggplot2')
library('plyr')
library("tidyverse")
options(bitmapType = 'cairo')

nc=nc_open(paste0(dir,'/Extremes.nc'))#for lonlat
nc_sv=nc_open(paste0(dir,'/Extremes_SV.nc'))#for lonlat
nc_obs=nc_open(paste0(dir,'/SeNorge.nc'))#for lonlat
nc_era=nc_open(paste0(dir,'/Extremes/ERA5.nc'))#for lonlat

Extremes_WC=ncvar_get(nc)
Extremes_SV=ncvar_get(nc_sv)
Extremes_obs=ncvar_get(nc_obs)
Extremes_era=ncvar_get(nc_era)
# Extremes_SV
dim(Extremes_WC) # 25 4 35 Ensemble Leadtime Year 

#There is an error that the dimnames do not get saved from Xarray to_netcdf. Set the dimnames here 
dimnames(Extremes_WC) = list(as.character(0:24),as.character(2:5),as.character(1981:2015))
dimnames(Extremes_SV) = list(as.character(0:24),as.character(2:5),as.character(1981:2015))
dimnames(Extremes_obs) = list(as.character(1957:2018))

predictand=as.vector(Extremes_obs[as.character(1981:2015)]) #First member, first leadtime that we use in this study
predictor=apply(Extremes_WC,MARGIN = c(2,3),FUN=mean) #predictor['2','1987']
predictor_all=apply(Extremes_WC,MARGIN = c(3),FUN=mean) #predictor['2','1987']

#Standarized anomaly
calc_anomaly <- function(variable) {
  (variable-mean(variable))/sd(variable)
}
predictor_anomaly=apply(predictor,MARGIN = 1 , FUN=calc_anomaly)
predictor_anomaly_all <- calc_anomaly(predictor_all)
predictand_anomaly=calc_anomaly(predictand)

