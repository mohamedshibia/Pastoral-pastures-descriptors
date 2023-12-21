#################################################################
##Passport Location
##GPS collar data import, file ordering and HHID, collar ID, time conversation
##Mohamed Shibia, mohamedshibia@gmail.com; s6moshib@uni-trier.de 
##Trier University, Germany
##################################################################
setwd("E:/Fern153/GPScollarData/cattledata")
#setwd("E:/GPScollarData/cattledata")#set working directly, CSV files folder
list.files() # name e.g. HHID_034_Collar_E037_Start_8-21-2011_End_3-24-2012
list_files <- list.files()
CattleTracks <- data.frame()

for(i in 1:length(list_files)){
  print(list_files[i])
  data <- read.csv(list_files[i]) 
  CattleTracks <- rbind(CattleTracks,data)
}

setwd("E:/Fern153/GPScollarData/outputs")
saveRDS(CattleTracks, "cattleGPSfiles.rds")

head(CattleTracks)# problem these data are not attributed to HHID, collarID
#write.csv(CattleTracks, file = "Cattle_GPSdata_All")

##write cattleTracks to R data frame data type
##Data carry individal confidential information, data available on request


#
setwd("E:/GPScollarData/outputs")
# 
# now <- proc.time()
# test2<-readRDS("cattlefile.rds")
# proc.time()-now

###identity a group and animal ID
setwd("E:/Fern153/GPScollarData/cattledata")
start <- proc.time()
newdata<-data.frame()
for(i in 1:length(list_files)){
  #print(list_files[i])
  print(sprintf("group: %s animal: %s",unlist(strsplit(list_files[i],"_"))[2],unlist(strsplit(list_files[i],"_"))[4]))
  data <- read.csv(list_files[i])
  groupcol<-rep(unlist(strsplit(list_files[i],"_"))[2],length(data[,1]))
  idcol<-rep(unlist(strsplit(list_files[i],"_"))[4],length(data[,1]))
  data<-cbind(data,groupcol,idcol)
  newdata <- rbind(newdata,data)
}
proc.time()-start

### confirmation testing using sub set by idcol
head(newdata)
a <- newdata[newdata$idcol=="E037",]
d1 <- newdata

saveRDS(newdata,file="cattleGPSfiles")
d1<-readRDS("cattleGPSfiles")#to read a RDS file in R
head(d1)
k <- d1[d1$idcol=="E037",] # raw data, groupid, idcol operational
## note we have more than 6 month trajectory per colar ID per year, 
#######################################################################################extracting elevation data at herd locations,calculating speed. Optional
###############################################################################
###############################################################################

require(rgdal)
require(raster)
#d1$Longitude[1],d1$Latitude[1]

#writing and calling a function
c(d1$Longitude[1],d1$Latitude[1])
getutm<-function(lon,lat){
  xy <- cbind(lon, lat)
  # Convert it to UTM coordinates (in units of meters)
  return(project(xy, "+proj=utm +zone=37N ellps=WGS84"))
}
getdist<-function(point1,point2){
  dv<-(c(point1[1]-point2[1],point1[2]-point2[2],point1[3]-point2[3]))
  return(sqrt(dv[1]^2+dv[2]^2+dv[3]^2))
}

#pass a matrix into the function
getutm(d1$Longitude[1],d1$Latitude[1])

#d1$Longitude<-as.numeric(d1$Longitude)
#d1$Elevation<-as.numeric(d1$Elevation)
#sum(is.na(d1$Longitude))
p1<-c(getutm(d1$Longitude[1],d1$Latitude[1]),d1$Elevation[1])
p2<-c(getutm(d1$Longitude[2],d1$Latitude[2]),d1$Elevation[2])
getdist(p1,p2)

d1$dist<-rep(0,length(d1[,1]))
for(i in 1:length(d1[,1])){
  p1<-c(getutm(d1$Longitude[i],d1$Latitude[i]),d1$Elevation[i])
  p2<-c(getutm(d1$Longitude[i+1],d1$Latitude[i+1]),d1$Elevation[i+1])
  d1$dist[i+1]<-getdist(p1,p2)
}

hist(d1$Elevation)
plot(d1$Elevation[1:1000],type="l")
# install.packages("signal")
# require(signal)
# 
# plot(d1$Elevation[1:1000],type="l")
# plot(sgolayfilt(d1$Elevation[1:1000],3,7,0,1),type="l")
dem<-raster("D:\\GPScollarData\\dem_msk_resampled") 
#dem<-raster("D:\\PAPER1\\A_DEM_mask.tif")
extent(dem)
crs(dem)
plot(dem)
extract(dem,cbind(500000,500000))

start <- proc.time() 
d1$Ele_DEM<-rep(0,length(d1[,1]))
for(i in 1:8322400){
  p<-getutm(d1$Longitude[i],d1$Latitude[i])
  d1$Ele_DEM[i]<-extract(dem,cbind(p[1],p[2]))
}
# g2=extract(dem,g1) #Chong 
# result=cbind(g1,g2)

proc.time()-start 
# for(i in 1:100){
#   p<-getutm(d1$Longitude[i],d1$Latitude[i])
#   d1$Ele_DEM[i]<-extract(dem,cbind(p[1],p[2]))
# }

## extract(dem,cbind(500000,500000))
## getutm(d1$Longitude[1],d1$Latitude[1])
## getutm(d1$Longitude[10],d1$Latitude[10])
## d1[10,]
saveRDS(d1, file = "d1_alldata_dem_elev.rds")

####
p1<-c(getutm(d1$Longitude[1],d1$Latitude[1]),d1$Ele_DEM[1])
p2<-c(getutm(d1$Longitude[2],d1$Latitude[2]),d1$Ele_DEM[2])
getdist(p1,p2)
d1$dist2<-rep(0,length(d1[,1]))
for(i in 1:length(d1[,1])){
  p1<-c(getutm(d1$Longitude[i],d1$Latitude[i]),d1$Ele_DEM[i])
  p2<-c(getutm(d1$Longitude[i+1],d1$Latitude[i+1]),d1$Ele_DEM[i+1])
  d1$dist2[i+1]<-getdist(p1,p2)
}
saveRDS(d1, file = "d1_alldata_dem_elev.rds") # OPTIONAL
######################################################################################################
####time difference calculation, speed calculation - REQUIRED 
#######################################################################################################
d1$UTC_Time
paste(d1$UTC_Date[1],d1$UTC_Time[1],sep=" ")
date1<-strptime(paste(d1$UTC_Date[1],d1$UTC_Time[1],sep=" "),format="%m/%d/%Y %H:%M:%S")
date2<-strptime(paste(d1$UTC_Date[2],d1$UTC_Time[2],sep=" "),format="%m/%d/%Y %H:%M:%S")

dt<-difftime(date2,date1,units="secs")
as.numeric(dt)

d1$deltat_seconds<-rep(0,length(d1[,1]))
d1$speed<-rep(0,length(d1[,1]))
for (i in 1:52352){
  date1<-strptime(paste(d1$UTC_Date[i],d1$UTC_Time[i],sep=" "),format="%m/%d/%Y %H:%M:%S")
  date2<-strptime(paste(d1$UTC_Date[i+1],d1$UTC_Time[i+1],sep=" "),format="%m/%d/%Y %H:%M:%S")
  dt<-difftime(date2,date1,units="secs")
  d1$deltat_seconds[i+1]<-as.numeric(dt)
  d1$speed<-d1$dist2/as.numeric(dt)
}

# ?difftime
# 
# as.difftime(c(as.character(d1$UTC_Time[1]),as.character(d1$UTC_Time[2])))

##############################################################################################################
####thresholding etc for data curation and confirmation test- OPTIONAL 
#############################################################################################################

par(mfrow=c(2,1))
plot(d1$Elevation[1:7200],type="l")
plot(d1$Ele_DEM[1:5000],type="l")
600*5
d1$date<-rep(as.POSIXct(time),length(d1[,1]))
for(i in 1:length(d1[,1])){
  time<-strptime(paste(d1$UTC_Date[i],d1$UTC_Time[i],sep=" "),format="%m/%d/%Y %H:%M:%S")
  d1$date[i]<-as.POSIXct(time)
  #d1$date[i]<-as.numeric(strftime(time,format="%H"))
  
}
par(mfrow=c(2,1))
start<-100
range<-40000
#plot(d1$date[1:range],d1$Elevation[1:range],type="l")
plot(d1$Ele_DEM[start:(start+range)],type="l")
plot(d1$date[start:(start+range)],d1$Ele_DEM[start:(start+range)],type="l")
plot(d1$Longitude,d1$Latitude,type="l")
hist(d1$dist2[1:500])
plot(d1$dist2,log="y")

# breaks<-c(seq(0,1000,by=200),seq(2000,100000,by=10000),max(na.exclude(d1$dist2))+1)
# hist(na.exclude(d1$dist2),breaks=500)
# max(na.exclude(d1$dist2))
thresh<-2000
subsetbelowthreshold<-d1[d1$dist2<thresh,]
length(d1[,1])-length(subsetbelowthreshold[,1])
length(subsetbelowthreshold[,1]) 

#write.csv(subsetbelowthreshold,file="IBLI_Cattle_GPS_COMPILED_Data_Feb 2015 to Aug 2015_sub")

