
# Open the netcdf files
# We should be in the UNSEEN-trends directory.
SEAS5_wc <- read.csv("Data/SEAS5_wc.csv") ## The UNSEEN ensemble
# SEAS5_sv_nc = nc_open('Data/Extremes_SV.nc')
ERA5 <- read.csv("Data/ERA5.csv") ## Reanalysis

head(SEAS5_wc)
SEAS5_wc$leadtime <- 10 - SEAS5_wc$leadtime
names(SEAS5_wc)[names(SEAS5_wc) == "leadtime"] <- "init-date"
names(SEAS5_wc)[names(SEAS5_wc) == "ensemble"] <- "member"

head(SEAS5_wc)
RV_20 <- stats::quantile(SEAS5_wc$precipitation, probs = 1-1/20) ## Select the quantile equal to return value 20-yr
SEAS5_wc_20yr <- SEAS5_wc[SEAS5_wc$precipitation >= RV_20,]


head(ERA5)
RV_10_era <- stats::quantile(ERA5$LSP, probs = 1-1/10) ## Select the quantile equal to return value 20-yr
ERA_10yr <- ERA5[ERA5$LSP >= RV_20_era,]


## And write the files
write.table(SEAS5_wc,file = '../../PhD/Work/UNSEEN-trends/Data/Malte/SEAS5.txt', row.names = FALSE)
write.table(SEAS5_wc_20yr[,1:3],file = '../../PhD/Work/UNSEEN-trends/Data/Malte/SEAS5_20yr.txt', row.names = FALSE)

write.table(ERA5,file = '../../PhD/Work/UNSEEN-trends/Data/Malte/ERA5.txt', row.names = FALSE)
write.table(ERA_10yr[,1],file = '../../PhD/Work/UNSEEN-trends/Data/Malte/ERA5_10yr.txt', row.names = FALSE)
