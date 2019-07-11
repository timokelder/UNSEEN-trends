library('ncdf4')
library('zoo')

filenames <- list.files(path = "//home/timok/timok/SALIENSEAS/SEAS5/ensex/Aggregated_months", pattern = glob2rx("seasonSON????.nc"),full.names = T)

years = format(seq(as.Date("1981/1/1"), as.Date("2015/1/1"), "years"),'%Y')
ensembles=4*25 #lead time * ensemble members

nc=nc_open(filenames[1])
longitude=ncvar_get(nc,'lon')
latitude=ncvar_get(nc,'lat')

lsp_array=ncvar_get(nc,'LSP') #lon,lat,time
lonlat <- as.matrix(expand.grid(longitude,latitude))
lsp_df01 <- data.frame(lonlat)
names(lsp_df01) <- c("lon","lat")
lsp_vec_long <- as.vector(lsp_array)
lsp_mat <- matrix(lsp_vec_long, nrow=length(longitude)*length(latitude), ncol=length(lsp_array[1,1,]))
lsp_df01[2+1:100]=lsp_mat

for (i in 2:length(years)){
  print(i)
  nc=nc_open(filenames[2])
  lsp_array=ncvar_get(nc,'LSP') #lon,lat,time
  lsp_vec_long <- as.vector(lsp_array)
  lsp_mat <- matrix(lsp_vec_long, nrow=length(longitude)*length(latitude), ncol=length(lsp_array[1,1,]))
  lsp_df01[2+(1+(i-1)*ensembles):(i*ensembles)]=lsp_mat
}

leadtimes=4
ld_list=rep(leadtimes:1+1,ensembles) #May-Aug forecasting Sep-Nov. Leadtime 1 is September, 2 = Aug.. May = 5. 
#May is sampled first, then june, july, aug and then next year. So the ld list is 5:2, repeated 35 times.

Extremes_ld2=lsp_df01[,3:length(lsp_df01[1,])][ld_list==2] #select only ld 2 from the ensembles, first two are latlon

Quantiles <- function(x) {quantile(x,probs=1-(1/200))}
a <-apply(Extremes_ld2,1,Quantiles)
lsp_df02 <- data.frame(lonlat)
lsp_df02[3]=a

write.table(lsp_df02,'//home/timok/timok/SALIENSEAS/SEAS5/ensex/statistics/multiday/Quantile_ld2/statistics_quantile200.txt' , row.names=FALSE, sep=",")


# create and write the netCDF file
# define dimensions
londim <- ncdim_def("lon","degrees_east",as.double(longitude)) 
latdim <- ncdim_def("lat","degrees_north",as.double(latitude)) 
# timedim <- ncdim_def("time",tunits3,as.double(time3))

# define variables
P_200_def <- ncvar_def("P_200","mm",list(londim,latdim))

# create netCDF file and put arrays
ncout <- nc_create('//home/timok/timok/SALIENSEAS/SEAS5/ensex/statistics/multiday/Quantile_ld2/Quantiles200.nc',P_200_def)

# put variables
P_200_array=array(lsp_df02[,3],dim= c(141,100))

ncvar_put(ncout,P_200_def,P_200_array)

# put additional attributes into dimension and data variables
ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")

# Get a summary of the created file:
ncout

# close the file, writing data to disk
nc_close(ncout)
