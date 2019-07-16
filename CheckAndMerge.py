import xarray as xr

dirname=r'//home/timok/timok/SALIENSEAS/SEAS5'

xr_year=xr.open_dataset(dirname+'/ensex/Extremes/Ens_lds1.nc')
for year in range(2,36): 
    xr_year_i=xr.open_dataset(dirname+'/ensex/Extremes/Ens_lds'+str(year)+'.nc')
    xr_year=xr.concat([xr_year,xr_year_i],'year')
    

xr_year=xr_year.assign_coords(year=range(1981,2016))


#Check the dataset: ensemble 1, Aug 2012 initialized =ld2
xr_1=xr.open_dataset('//home/timok/timok/SALIENSEAS/SEAS5/01/Arctic.SEAS5_sfc_01_20120801.nc')
xr1=xr_1.sel(time=slice('2012-9-01','2012-12-01'),lon=slice(4,7),lat=slice(58,63)).drop(['SSTK','CI']) #drop the unnecessary variables
#open the netcdf with 200-year return quantiles to mask regions
xr_q200=xr.open_dataset(dirname+'/ensex/statistics/multiday/Quantile_ld2/Quantiles200.nc')
#The region is selected based on where the climatology (200-year values) are greater than a user-defined threshold
climatology_threshold=90
#Select the regional averaged 3 day cumulative precipitation
Extreme_ens01_ld_02_2012=xr1['LSP'].diff('time').rolling(time=3).sum().where(xr_q200['P_200']>climatology_threshold).mean(dim=['lat','lon']).max(dim='time')*1000
xr_year.sel(ensemble=1,leadtime=2,year=2012).values()==Extreme_ens01_ld_02_2012.values #True :)

xr_year.to_netcdf('//home/timok/timok/SALIENSEAS/SEAS5/ensex/Extremes/Extremes.nc', mode='w')

