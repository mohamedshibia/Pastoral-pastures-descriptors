#https://trinkerrstuff.wordpress.com/2013/09/15/paste-paste0-and-sprintf-2/ 

#install.packages("lubridate")
library(adehabitatLT)
library(lubridate)
library(raster)
library(rgdal)
library(hab)
library(sp)


#all raw dataset for round2
#data<-readRDS("Round2.rds")

###sub herd selection
###############################################################################################################
############################################# hhid 730, name Bidu Wario Guyo###################################
###############################################################################################################
#pastoralist data
#setwd("D:/GPScollarData/cattledata")
#setwd("D:/shibia/herds") 
list.files()
list_files <- list.files()
length <- length(list_files)

##write cattleTracks to R data frame data type
###identiy a group and animal
start <- proc.time()
newdata<-data.frame()# pre-allocate some variable space for our output data
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

head(newdata)
p.730 <- newdata
saveRDS(newdata,file="p.730.allrounds.data.raw.rds") 

sum(is.na(newdata$Longitude))
sum(is.na(newdata$Latitude))
sum(is.na(newdata$Elevation))
sum(is.na(newdata$UTC_Date))
sum(is.na(newdata$UTC_Time))


plot(newdata$Longitude, newdata$Latitude, col=newdata$idcol) #pastoralist stayed in a permanent place throughout study period.

da<-strptime(paste(newdata$UTC_Date,newdata$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #mario 

#r1
r1start<-as.POSIXct("2011-08-01 00:00:00 CEST")
r1end<-as.POSIXct("2012-03-31 00:00:00 CEST")
#r2
r2start<-as.POSIXct("2012-09-01 00:00:00 CEST")
r2end<-as.POSIXct("2013-02-27 00:00:00 CEST")
#r3
r3start<-as.POSIXct("2013-09-01 00:00:00 CEST")
r3end<-as.POSIXct("2014-02-14 00:00:00 CEST")
#r4
r4start<-as.POSIXct("2014-02-15 00:00:00 CEST")
r4end<-as.POSIXct("2014-08-31 00:00:00 CEST")
#r5
r5start<-as.POSIXct("2014-09-01 00:00:00 CEST")
r5end<-as.POSIXct("2015-02-12 00:00:00 CEST")
#r6
r6start<-as.POSIXct("2015-02-13 00:00:00 CEST")#not 15th feb but 13th feb, and contains about 292 NAs
r6end<-as.POSIXct("2015-08-31 00:00:00 CEST")

da[2]<da[1]
classround<-function(date){
  if (is.na(date)) return (NA)
  if (r1start<date & r1end> date) return("r1")
  if (r2start<date & r2end> date) return("r2")
  if (r3start<date & r3end> date) return("r3")
  if (r4start<date & r4end> date) return("r4")
  if (r5start<date & r5end> date) return("r5")
  if (r6start<date & r6end> date) return("r6")
  else return (NA)
}

#da_class<-lapply(da,classround) #not working 
#foo to classify da and fill in a vector 
dclass<-vector()
for (i in 1:length(da)){
  print(i)
  dclass[i]<-classround(da[i])
}
daclass_f<-as.factor(dclass) #curse class as.factor 
summary(daclass_f) #summary 

summary(daclass_f[newdata$idcol=="E049"]) # a test 
for( i in 1:length(unique(newdata$idcol))){
  id<-unique(newdata$idcol)[i]
  print(id)
  print(summary(daclass_f[newdata$idcol==id]))
}


#######################################################################################################################
######################################################################################################################
#######################################################################################################################
######################################################################################################################
#select an ID of GPS collar that at least recorded sufficient information for an animal group in every season.
#create a new variable called Rounds in the data frame to identify relocation points based on the UTC_Date  

newdata$Rounds <- daclass_f
table(newdata$idcol, newdata$Rounds)



#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E051' & newdata$Rounds=='r1'),] #60195 obs.
newdata.r2 <- newdata[which(newdata$idcol=='E051' & newdata$Rounds=='r2'),]#45523
newdata.r3 <- newdata[which(newdata$idcol=='e050' & newdata$Rounds=='r3'),]# 4115
newdata.r4 <- newdata[which(newdata$idcol=='e112' & newdata$Rounds=='r4'),]#49564
newdata.r5 <- newdata[which(newdata$idcol=='e112' & newdata$Rounds=='r5'),]#5144
newdata.r6 <- newdata[which(newdata$idcol=='e103' & newdata$Rounds=='r6'),]# 24727

hhid.730.herd.w <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.730.herd.w)
plot(hhid.730.herd.w$Longitude, hhid.730.herd.w$Latitude, col=hhid.730.herd.w$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  
hhid.730.herd.w.sub<- subset(hhid.730.herd.w, hhid.730.herd.w$Longitude>38.2 & hhid.730.herd.w$Longitude<38.5 & hhid.730.herd.w$Latitude>4.5 & hhid.730.herd.w$Latitude<5.0)
plot(hhid.730.herd.w.sub$Longitude, hhid.730.herd.w.sub$Latitude, col=hhid.730.herd.w.sub$Rounds)

#create a uniquely new ID for these collars.  
hhid.730.herd.w.sub$idcol_new <- as.factor("c001")

#ensure files are of equal lengths using key variables  
length(hhid.730.herd.w.sub$idcol_new) 
length(hhid.730.herd.w.sub$Longitude)
length(hhid.730.herd.w.sub$Latitude) 
length(hhid.730.herd.w.sub$UTC_Date) 
length(hhid.730.herd.w.sub$UTC_Time) 

saveRDS(hhid.730.herd.w.sub, file = "c001.hhid.730.wara.subherd.rds")
newdata <- readRDS("c001.hhid.730.wara.subherd.rds")
plot(newdata$Longitude, newdata$Latitude, col=newdata$Rounds)
summary(newdata$idcol_new)

##a file clean up
newdata.sub<- subset(newdata, newdata$Longitude>38.34 & newdata$Longitude<38.42 & newdata$Latitude>4.5 & newdata$Latitude<5.0)
plot(newdata.sub$Longitude, newdata.sub$Latitude, col=newdata.sub$Rounds)
summary(newdata.sub$idcol_new)
saveRDS(newdata.sub, file = "c001.hhid.730.wara.subherd.rds")

#################################################################################################################
#################################################################################################################

###hhid 774, name Dika Garbicha Jatani

################################################################################################################
################################################################################################################
#import csv files 

list.files()
list_files <- list.files()
length <- length(list_files)

##write cattleTracks to R data frame data type
###identiy a group and animal
start <- proc.time()
newdata<-data.frame()# pre-allocate some variable space for our output data
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

head(newdata)
#p.730 <- newdata
saveRDS(newdata,file="p.774.allrounds.data.raw.rds") 
newdata <- readRDS("p.774.allrounds.data.raw.rds")


sum(is.na(newdata$Longitude))
sum(is.na(newdata$Latitude))
sum(is.na(newdata$Elevation))
sum(is.na(newdata$UTC_Date))
sum(is.na(newdata$UTC_Time))


plot(newdata$Longitude, newdata$Latitude, col=newdata$idcol) #pastoralist stayed in a permanent place throughout study period.

da<-strptime(paste(newdata$UTC_Date,newdata$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #mario 

#r1
r1start<-as.POSIXct("2011-08-01 00:00:00 CEST")
r1end<-as.POSIXct("2012-03-31 00:00:00 CEST")
#r2
r2start<-as.POSIXct("2012-09-01 00:00:00 CEST")
r2end<-as.POSIXct("2013-02-27 00:00:00 CEST")
#r3
r3start<-as.POSIXct("2013-09-01 00:00:00 CEST")
r3end<-as.POSIXct("2014-02-14 00:00:00 CEST")
#r4
r4start<-as.POSIXct("2014-02-15 00:00:00 CEST")
r4end<-as.POSIXct("2014-08-31 00:00:00 CEST")
#r5
r5start<-as.POSIXct("2014-09-01 00:00:00 CEST")
r5end<-as.POSIXct("2015-02-12 00:00:00 CEST")
#r6
r6start<-as.POSIXct("2015-02-13 00:00:00 CEST")#not 15th feb but 13th feb, and contains about 292 NAs
r6end<-as.POSIXct("2015-08-31 00:00:00 CEST")

da[2]<da[1]
classround<-function(date){
  if (is.na(date)) return (NA)
  if (r1start<date & r1end> date) return("r1")
  if (r2start<date & r2end> date) return("r2")
  if (r3start<date & r3end> date) return("r3")
  if (r4start<date & r4end> date) return("r4")
  if (r5start<date & r5end> date) return("r5")
  if (r6start<date & r6end> date) return("r6")
  else return (NA)
}

#da_class<-lapply(da,classround) #not working 
#foo to classify da and fill in a vector 
dclass<-vector()
for (i in 1:length(da)){
  print(i)
  dclass[i]<-classround(da[i])
}
daclass_f<-as.factor(dclass) #curse class as.factor 
summary(daclass_f) #summary 

#summary(daclass_f[newdata$idcol=="E049"]) # a test 
for( i in 1:length(unique(newdata$idcol))){
  id<-unique(newdata$idcol)[i]
  print(id)
  print(summary(daclass_f[newdata$idcol==id]))
}


#######################################################################################################################
######################################################################################################################
#######################################################################################################################
######################################################################################################################
#select an ID of GPS collar that at least recorded sufficient information for an animal group in every season.
#create a new variable called Rounds in the data frame to identify relocation points based on the UTC_Date  

newdata$Rounds <- daclass_f
head(newdata)
table(newdata$idcol, newdata$Rounds)

#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E054' & newdata$Rounds=='r1'),] #46903 obs.
newdata.r2 <- newdata[which(newdata$idcol=='E054' & newdata$Rounds=='r2'),]#45523
newdata.r3 <- newdata[which(newdata$idcol=='e052' & newdata$Rounds=='r3'),]# 43837
newdata.r4 <- newdata[which(newdata$idcol=='e052' & newdata$Rounds=='r4'),]#49608
newdata.r5 <- newdata[which(newdata$idcol=='e108' & newdata$Rounds=='r5'),]#40323
newdata.r6 <- newdata[which(newdata$idcol=='e052' & newdata$Rounds=='r6'),]# 25239

hhid.774.herd.w <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.774.herd.w)
plot(hhid.774.herd.w$Longitude, hhid.774.herd.w$Latitude, col=hhid.774.herd.w$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  
hhid.774.herd.w.sub<- subset(hhid.774.herd.w, hhid.774.herd.w$Longitude>38.36 & hhid.774.herd.w$Longitude<38.46 & hhid.774.herd.w$Latitude>4.6 & hhid.774.herd.w$Latitude<4.86)
plot(hhid.774.herd.w.sub$Longitude, hhid.774.herd.w.sub$Latitude, col=hhid.774.herd.w.sub$Rounds)

#create a uniquely new ID for these collars.  
hhid.774.herd.w.sub$idcol_new <- as.factor("c002")

#ensure equal length for key variables  
length(hhid.774.herd.w.sub$idcol_new) 
length(hhid.774.herd.w.sub$Longitude)
length(hhid.774.herd.w.sub$Latitude) 
length(hhid.774.herd.w.sub$UTC_Date) 
length(hhid.774.herd.w.sub$UTC_Time) 



saveRDS(hhid.774.herd.w.sub, file = "c002.hhid.774.wara.subherd.rds")

newdata <- readRDS("c002.hhid.774.wara.subherd.rds")
plot(newdata$Longitude, newdata$Latitude, col=newdata$Rounds)
summary(newdata$idcol_new)



#################################################################################################################
#################################################################################################################

###hhid 807, name Gufuu Galgallo duubaa 

################################################################################################################
################################################################################################################
#import csv files 

list.files()
list_files <- list.files()
length <- length(list_files)

##write cattleTracks to R data frame data type
###identiy a group and animal
start <- proc.time()
newdata<-data.frame()# pre-allocate some variable space for our output data
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

head(newdata)
#p.730 <- newdata
saveRDS(newdata,file="p.807.allrounds.data.raw.rds") 
newdata <- readRDS("p.807.allrounds.data.raw.rds")

sum(is.na(newdata$Longitude))
sum(is.na(newdata$Latitude))
sum(is.na(newdata$Elevation))
sum(is.na(newdata$UTC_Date))
sum(is.na(newdata$UTC_Time))


plot(newdata$Longitude, newdata$Latitude, col=newdata$idcol) #pastoralist stayed in a permanent and temporary place.

da<-strptime(paste(newdata$UTC_Date,newdata$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #mario 

#r1
r1start<-as.POSIXct("2011-08-01 00:00:00 CEST")
r1end<-as.POSIXct("2012-03-31 00:00:00 CEST")
#r2
r2start<-as.POSIXct("2012-09-01 00:00:00 CEST")
r2end<-as.POSIXct("2013-02-27 00:00:00 CEST")
#r3
r3start<-as.POSIXct("2013-09-01 00:00:00 CEST")
r3end<-as.POSIXct("2014-02-14 00:00:00 CEST")
#r4
r4start<-as.POSIXct("2014-02-15 00:00:00 CEST")
r4end<-as.POSIXct("2014-08-31 00:00:00 CEST")
#r5
r5start<-as.POSIXct("2014-09-01 00:00:00 CEST")
r5end<-as.POSIXct("2015-02-12 00:00:00 CEST")
#r6
r6start<-as.POSIXct("2015-02-13 00:00:00 CEST")#not 15th feb but 13th feb, and contains about 292 NAs
r6end<-as.POSIXct("2015-09-30 00:00:00 CEST") #changed to sep 30 from aug 31

da[2]<da[1]
classround<-function(date){
  if (is.na(date)) return (NA)
  if (r1start<date & r1end> date) return("r1")
  if (r2start<date & r2end> date) return("r2")
  if (r3start<date & r3end> date) return("r3")
  if (r4start<date & r4end> date) return("r4")
  if (r5start<date & r5end> date) return("r5")
  if (r6start<date & r6end> date) return("r6")
  else return (NA)
}

#da_class<-lapply(da,classround) #not working 
#foo to classify da and fill in a vector 
dclass<-vector()
for (i in 1:length(da)){
  print(i)
  dclass[i]<-classround(da[i])
}
daclass_f<-as.factor(dclass) #curse class as.factor 
summary(daclass_f) #summary 

#summary(daclass_f[newdata$idcol=="E049"]) # a test 
for( i in 1:length(unique(newdata$idcol))){
  id<-unique(newdata$idcol)[i]
  print(id)
  print(summary(daclass_f[newdata$idcol==id]))
}

#######################################################################################################################
######################################################################################################################
#select an ID of GPS collar that at least recorded sufficient information for an animal group in every season.
#create a new variable called Rounds in the data frame to identify relocation points based on the UTC_Date  

newdata$Rounds <- daclass_f
head(newdata)
table(newdata$idcol, newdata$Rounds)

#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E060' & newdata$Rounds=='r1'),] #56604 obs.
newdata.r2 <- newdata[which(newdata$idcol=='E059' & newdata$Rounds=='r2'),]#45577
newdata.r3 <- newdata[which(newdata$idcol=='e058' & newdata$Rounds=='r3'),]# 43837
newdata.r4 <- newdata[which(newdata$idcol=='e104' & newdata$Rounds=='r4'),]#49638
newdata.r5 <- newdata[which(newdata$idcol=='e104' & newdata$Rounds=='r5'),]#41271
newdata.r6 <- newdata[which(newdata$idcol=='e104' & newdata$Rounds=='r6'),]# 35069


hhid.807.herd.f <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.807.herd.f)
plot(hhid.807.herd.f$Longitude, hhid.807.herd.f$Latitude, col=hhid.807.herd.f$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  
hhid.807.herd.f.sub<- subset(hhid.807.herd.f, hhid.807.herd.f$Longitude>38.2 & hhid.807.herd.f$Longitude<39.5 & hhid.807.herd.f$Latitude>4.0 & hhid.807.herd.f$Latitude<4.86)
plot(hhid.807.herd.f.sub$Longitude, hhid.807.herd.f.sub$Latitude, col=hhid.807.herd.f.sub$Rounds)

#create a uniquely new ID for these collars.  
hhid.807.herd.f.sub$idcol_new <- as.factor("c003")

#ensure equal length for key variables  
length(hhid.807.herd.f.sub$idcol_new) 
length(hhid.807.herd.f.sub$Longitude)
length(hhid.807.herd.f.sub$Latitude) 
length(hhid.807.herd.f.sub$UTC_Date) 
length(hhid.807.herd.f.sub$UTC_Time) 



saveRDS(hhid.807.herd.f.sub, file = "c003.hhid.807.fora.subherd.rds")
newdata <- readRDS("c003.hhid.807.fora.subherd.rds")
plot(newdata$Longitude, newdata$Latitude, col=newdata$Rounds)
summary(newdata$idcol_new)
 

##create a new data for wara subherd 


newdata.r1 <-newdata[which(newdata$idcol=='E059' & newdata$Rounds=='r1'),] #56604 obs.
newdata.r2 <- newdata[which(newdata$idcol=='E060' & newdata$Rounds=='r2'),]#45577
#newdata.r3 <- newdata[which(newdata$idcol=='e058' & newdata$Rounds=='r3'),]# 43837
newdata.r4 <- newdata[which(newdata$idcol=='e116' & newdata$Rounds=='r4'),]#49638
newdata.r5 <- newdata[which(newdata$idcol=='e116' & newdata$Rounds=='r5'),]#41271
newdata.r6 <- newdata[which(newdata$idcol=='e116' & newdata$Rounds=='r6'),]# 35069


hhid.807.herd.w <- rbind(newdata.r1, newdata.r2, newdata.r4, newdata.r5, newdata.r6)
head(hhid.807.herd.w)
plot(hhid.807.herd.w$Longitude, hhid.807.herd.w$Latitude, col=hhid.807.herd.w$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  
# hhid.807.herd.f.sub<- subset(hhid.807.herd.f, hhid.807.herd.f$Longitude>38.2 & hhid.807.herd.f$Longitude<39.5 & hhid.807.herd.f$Latitude>4.0 & hhid.807.herd.f$Latitude<4.86)
# plot(hhid.807.herd.f.sub$Longitude, hhid.807.herd.f.sub$Latitude, col=hhid.807.herd.f.sub$Rounds)

#create a uniquely new ID for these collars.  
hhid.807.herd.w$idcol_new <- as.factor("c003.w")

#ensure equal length for key variables  
length(hhid.807.herd.w$idcol_new) 
length(hhid.807.herd.w$Longitude)
length(hhid.807.herd.w$Latitude) 
length(hhid.807.herd.w$UTC_Date) 
length(hhid.807.herd.w$UTC_Time) 



saveRDS(hhid.807.herd.w, file = "c003.w.hhid.807.wora.subherd.rds") #wara sub herd 




#################################################################################################################
#################################################################################################################

###hhid 139, name Nura Guyo Ilme  

################################################################################################################
################################################################################################################
#import csv files 

list.files()
list_files <- list.files()
length <- length(list_files)

##write cattleTracks to R data frame data type
###identiy a group and animal
start <- proc.time()
newdata<-data.frame()# pre-allocate some variable space for our output data
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

head(newdata)
#p.730 <- newdata
saveRDS(newdata,file="p.139.allrounds.data.raw.rds") 

sum(is.na(newdata$Longitude))
sum(is.na(newdata$Latitude))
sum(is.na(newdata$Elevation))
sum(is.na(newdata$UTC_Date))
sum(is.na(newdata$UTC_Time))


plot(newdata$Longitude, newdata$Latitude, col=newdata$idcol) #pastoralist had both temporary and permanent place since the inception of study.

da<-strptime(paste(newdata$UTC_Date,newdata$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #mario 

#r1
r1start<-as.POSIXct("2011-08-01 00:00:00 CEST")
r1end<-as.POSIXct("2012-03-31 00:00:00 CEST")
#r2
r2start<-as.POSIXct("2012-09-01 00:00:00 CEST")
r2end<-as.POSIXct("2013-02-27 00:00:00 CEST")
#r3
r3start<-as.POSIXct("2013-09-01 00:00:00 CEST")
r3end<-as.POSIXct("2014-02-14 00:00:00 CEST")
#r4
r4start<-as.POSIXct("2014-02-15 00:00:00 CEST")
r4end<-as.POSIXct("2014-08-31 00:00:00 CEST")
#r5
r5start<-as.POSIXct("2014-09-01 00:00:00 CEST")
r5end<-as.POSIXct("2015-02-12 00:00:00 CEST")
#r6
r6start<-as.POSIXct("2015-02-13 00:00:00 CEST")#not 15th feb but 13th feb, and contains about 292 NAs
r6end<-as.POSIXct("2015-09-30 00:00:00 CEST") #changed to sep 30 from aug 31

da[2]<da[1]
classround<-function(date){
  if (is.na(date)) return (NA)
  if (r1start<date & r1end> date) return("r1")
  if (r2start<date & r2end> date) return("r2")
  if (r3start<date & r3end> date) return("r3")
  if (r4start<date & r4end> date) return("r4")
  if (r5start<date & r5end> date) return("r5")
  if (r6start<date & r6end> date) return("r6")
  else return (NA)
}

#da_class<-lapply(da,classround) #not working 
#foo to classify da and fill in a vector 
dclass<-vector()
for (i in 1:length(da)){
  print(i)
  dclass[i]<-classround(da[i])
}
daclass_f<-as.factor(dclass) #curse class as.factor 
summary(daclass_f) #summary 

#summary(daclass_f[newdata$idcol=="E049"]) # a test 
for( i in 1:length(unique(newdata$idcol))){
  id<-unique(newdata$idcol)[i]
  print(id)
  print(summary(daclass_f[newdata$idcol==id]))
}

#######################################################################################################################
######################################################################################################################
#select an ID of GPS collar that at least recorded sufficient information for an animal group in every season.
#create a new variable called Rounds in the data frame to identify relocation points based on the UTC_Date  

newdata$Rounds <- daclass_f
head(newdata)
table(newdata$idcol, newdata$Rounds)

#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E025' & newdata$Rounds=='r1'),] #57801 obs.
newdata.r2 <- newdata[which(newdata$idcol=='E031' & newdata$Rounds=='r2'),]#44460
newdata.r3 <- newdata[which(newdata$idcol=='e026' & newdata$Rounds=='r3'),]# 38184
newdata.r4 <- newdata[which(newdata$idcol=='e031' & newdata$Rounds=='r4'),]#47133
newdata.r5 <- newdata[which(newdata$idcol=='e031' & newdata$Rounds=='r5'),]#41142
newdata.r6 <- newdata[which(newdata$idcol=='e031' & newdata$Rounds=='r6'),]# 17623


hhid.139.herd.f <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.139.herd.f)
plot(hhid.139.herd.f$Longitude, hhid.139.herd.f$Latitude, col=hhid.139.herd.f$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  
hhid.139.herd.f.sub<- subset(hhid.139.herd.f, hhid.139.herd.f$Longitude>38.5 & hhid.139.herd.f$Longitude<39.0 & hhid.139.herd.f$Latitude>3.8 & hhid.139.herd.f$Latitude<4.3)
plot(hhid.139.herd.f.sub$Longitude, hhid.139.herd.f.sub$Latitude, col=hhid.139.herd.f.sub$Rounds)

#create a uniquely new ID for these collars.  
hhid.139.herd.f.sub$idcol_new <- as.factor("c004")

#ensure equal length for key variables  
length(hhid.139.herd.f.sub$idcol_new) 
length(hhid.139.herd.f.sub$Longitude)
length(hhid.139.herd.f.sub$Latitude) 
length(hhid.139.herd.f.sub$UTC_Date) 
length(hhid.139.herd.f.sub$UTC_Time) 



saveRDS(hhid.139.herd.f.sub, file = "c004.hhid.139.fora.subherd.rds")

newdata <- readRDS("c004.hhid.139.fora.subherd.rds")
plot(newdata$Longitude, newdata$Latitude, col=newdata$Rounds)
summary(newdata$idcol_new)



############### wara sub herd for hhig 139
#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E026' & newdata$Rounds=='r1'),] #23173 obs.
newdata.r2 <- newdata[which(newdata$idcol=='E026' & newdata$Rounds=='r2'),]#44518
newdata.r3 <- newdata[which(newdata$idcol=='e065' & newdata$Rounds=='r3'),]# 38118
newdata.r4 <- newdata[which(newdata$idcol=='e026' & newdata$Rounds=='r4'),]#8917
newdata.r5 <- newdata[which(newdata$idcol=='e026' & newdata$Rounds=='r5'),]#40706
newdata.r6 <- newdata[which(newdata$idcol=='e026' & newdata$Rounds=='r6'),]# 15140


hhid.139.herd.w <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.139.herd.w)
plot(hhid.139.herd.w$Longitude, hhid.139.herd.w$Latitude, col=hhid.139.herd.w$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  
hhid.139.herd.w.sub<- subset(hhid.139.herd.w, hhid.139.herd.w$Longitude>38.55 & hhid.139.herd.w$Longitude<38.85 & hhid.139.herd.w$Latitude>3.9 & hhid.139.herd.w$Latitude<4.3)
plot(hhid.139.herd.w.sub$Longitude, hhid.139.herd.w.sub$Latitude, col=hhid.139.herd.w.sub$Rounds)

#create a uniquely new ID for these collars.  
hhid.139.herd.w.sub$idcol_new <- as.factor("c004.w")

#ensure equal length for key variables  
length(hhid.139.herd.w.sub$idcol_new) 
length(hhid.139.herd.w.sub$Longitude)
length(hhid.139.herd.w.sub$Latitude) 
length(hhid.139.herd.w.sub$UTC_Date) 
length(hhid.139.herd.w.sub$UTC_Time) 



saveRDS(hhid.139.herd.w.sub, file = "c004.w.hhid.139.wara.subherd.rds") #wara sub herd 

newdata <- readRDS("c004.w.hhid.139.wara.subherd.rds")
plot(newdata$Longitude, newdata$Latitude, col=newdata$Rounds)
summary(newdata$idcol_new)


#################################################################################################################
#################################################################################################################

###hhid 155, name Roba Kulu waqo  

################################################################################################################
################################################################################################################
#import csv files 

list.files()
list_files <- list.files()
length <- length(list_files)

##write cattleTracks to R data frame data type
###identiy a group and animal
start <- proc.time()
newdata<-data.frame()# pre-allocate some variable space for our output data
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

head(newdata)
#p.730 <- newdata
saveRDS(newdata,file="p.155.allrounds.data.raw.rds") 
newdata <- readRDS("p.155.allrounds.data.raw.rds")

sum(is.na(newdata$Longitude)) #a total of 0 nas
sum(is.na(newdata$Latitude)) #a total of 0 nas 
sum(is.na(newdata$Elevation)) #a total of 0 nas
sum(is.na(newdata$UTC_Date))
sum(is.na(newdata$UTC_Time))

#select complete data without nas
# n.data <- readRDS("p.155.allrounds.data.raw.rds") 
# n.data.2 <- n.data[!is.na(n.data$Longitude),]
# head(n.data.2)
# sum(is.na(n.data.2$Longitude)) #a total of 0 nas
# sum(is.na(n.data.2$Latitude)) #a total of 0 nas 
# sum(is.na(n.data.2$Elevation)) 
# sum(is.na(n.data.2$UTC_Date))
# sum(is.na(n.data.2$UTC_Time))

#new_DF <- DF[is.na(DF$Var),] #an example



plot(newdata$Longitude, newdata$Latitude, col=newdata$idcol) #pastoralist had both temporary and permanent place since the inception of study.

da<-strptime(paste(newdata$UTC_Date,newdata$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #mario 

#r1
r1start<-as.POSIXct("2011-08-01 00:00:00 CEST")
r1end<-as.POSIXct("2012-03-31 00:00:00 CEST")
#r2
r2start<-as.POSIXct("2012-09-01 00:00:00 CEST")
r2end<-as.POSIXct("2013-02-27 00:00:00 CEST")
#r3
r3start<-as.POSIXct("2013-09-01 00:00:00 CEST")
r3end<-as.POSIXct("2014-02-14 00:00:00 CEST")
#r4
r4start<-as.POSIXct("2014-02-15 00:00:00 CEST")
r4end<-as.POSIXct("2014-08-31 00:00:00 CEST")
#r5
r5start<-as.POSIXct("2014-09-01 00:00:00 CEST")
r5end<-as.POSIXct("2015-02-12 00:00:00 CEST")
#r6
r6start<-as.POSIXct("2015-02-13 00:00:00 CEST")#not 15th feb but 13th feb, and contains about 292 NAs
r6end<-as.POSIXct("2015-09-30 00:00:00 CEST") #changed to sep 30 from aug 31

da[2]<da[1]
classround<-function(date){
  if (is.na(date)) return (NA)
  if (r1start<date & r1end> date) return("r1")
  if (r2start<date & r2end> date) return("r2")
  if (r3start<date & r3end> date) return("r3")
  if (r4start<date & r4end> date) return("r4")
  if (r5start<date & r5end> date) return("r5")
  if (r6start<date & r6end> date) return("r6")
  else return (NA)
}

#da_class<-lapply(da,classround) #not working 
#foo to classify da and fill in a vector 
dclass<-vector()
for (i in 1:length(da)){
  print(i)
  dclass[i]<-classround(da[i])
}
daclass_f<-as.factor(dclass) #curse class as.factor 
summary(daclass_f) #summary 

#summary(daclass_f[newdata$idcol=="E049"]) # a test 
for( i in 1:length(unique(newdata$idcol))){
  id<-unique(newdata$idcol)[i]
  print(id)
  print(summary(daclass_f[newdata$idcol==id]))
}

#######################################################################################################################
######################################################################################################################
#select an ID of GPS collar that at least recorded sufficient information for an animal group in every season.
#create a new variable called Rounds in the data frame to identify relocation points based on the UTC_Date  

newdata$Rounds <- daclass_f
head(newdata)
table(newdata$idcol, newdata$Rounds)

#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E033' & newdata$Rounds=='r1'),] #55812 obs.
newdata.r2 <- newdata[which(newdata$idcol=='E033' & newdata$Rounds=='r2'),]#44644
newdata.r3 <- newdata[which(newdata$idcol=='e033' & newdata$Rounds=='r3'),]# 38036
newdata.r4 <- newdata[which(newdata$idcol=='r067' & newdata$Rounds=='r4'),]#47184
newdata.r5 <- newdata[which(newdata$idcol=='e033' & newdata$Rounds=='r5'),]#13770
newdata.r6 <- newdata[which(newdata$idcol=='e033' & newdata$Rounds=='r6'),]# 17309


hhid.155.herd.f <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.155.herd.f)
plot(hhid.155.herd.f$Longitude, hhid.155.herd.f$Latitude, col=hhid.155.herd.f$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  
hhid.155.herd.f.sub<- subset(hhid.155.herd.f, hhid.155.herd.f$Longitude>38.55 & hhid.155.herd.f$Longitude<38.75 & hhid.155.herd.f$Latitude>3.7 & hhid.155.herd.f$Latitude<4.13)
plot(hhid.155.herd.f.sub$Longitude, hhid.155.herd.f.sub$Latitude, col=hhid.155.herd.f.sub$Rounds)

#create a uniquely new ID for these collars.  
hhid.155.herd.f.sub$idcol_new <- as.factor("c005")

#ensure equal length for key variables  
length(hhid.155.herd.f.sub$idcol_new) 
length(hhid.155.herd.f.sub$Longitude)
length(hhid.155.herd.f.sub$Latitude) 
length(hhid.155.herd.f.sub$UTC_Date) 
length(hhid.155.herd.f.sub$UTC_Time) 



saveRDS(hhid.155.herd.f.sub, file = "c005.hhid.155.fora.subherd.rds")

newdata <- readRDS("c005.hhid.155.fora.subherd.rds")
plot(newdata$Longitude, newdata$Latitude, col=newdata$Rounds)
summary(newdata$idcol_new)


#################################################################################################################
#################################################################################################################

###hhid 184, name Diida Duba Sora 

################################################################################################################
################################################################################################################
#import csv files 

list.files()
list_files <- list.files()
length <- length(list_files)

##write cattleTracks to R data frame data type
###identiy a group and animal
start <- proc.time()
newdata<-data.frame()# pre-allocate some variable space for our output data
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

head(newdata)
#p.730 <- newdata
saveRDS(newdata,file="p.184.allrounds.data.raw.rds") 

sum(is.na(newdata$Longitude)) #a total of 0 nas
sum(is.na(newdata$Latitude)) #a total of 0 nas 
sum(is.na(newdata$Elevation)) #a total of 0 nas
sum(is.na(newdata$UTC_Date))
sum(is.na(newdata$UTC_Time))

#choose to run these lines only to select complete data without nas
# n.data <- readRDS("p.155.allrounds.data.raw.rds") 
# n.data.2 <- n.data[!is.na(n.data$Longitude),]
# head(n.data.2)
# sum(is.na(n.data.2$Longitude)) #a total of 0 nas
# sum(is.na(n.data.2$Latitude)) #a total of 0 nas 
# sum(is.na(n.data.2$Elevation)) 
# sum(is.na(n.data.2$UTC_Date))
# sum(is.na(n.data.2$UTC_Time))

#new_DF <- DF[is.na(DF$Var),] #an example



plot(newdata$Longitude, newdata$Latitude, col=newdata$idcol) #pastoralist had both temporary and permanent place since the inception of study.

da<-strptime(paste(newdata$UTC_Date,newdata$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #mario 

#r1
r1start<-as.POSIXct("2011-08-01 00:00:00 CEST")
r1end<-as.POSIXct("2012-03-31 00:00:00 CEST")
#r2
r2start<-as.POSIXct("2012-09-01 00:00:00 CEST")
r2end<-as.POSIXct("2013-02-27 00:00:00 CEST")
#r3
r3start<-as.POSIXct("2013-09-01 00:00:00 CEST")
r3end<-as.POSIXct("2014-02-14 00:00:00 CEST")
#r4
r4start<-as.POSIXct("2014-02-15 00:00:00 CEST")
r4end<-as.POSIXct("2014-08-31 00:00:00 CEST")
#r5
r5start<-as.POSIXct("2014-09-01 00:00:00 CEST")
r5end<-as.POSIXct("2015-02-12 00:00:00 CEST")
#r6
r6start<-as.POSIXct("2015-02-13 00:00:00 CEST")#not 15th feb but 13th feb, and contains about 292 NAs
r6end<-as.POSIXct("2015-09-30 00:00:00 CEST") #changed to sep 30 from aug 31

da[2]<da[1]
classround<-function(date){
  if (is.na(date)) return (NA)
  if (r1start<date & r1end> date) return("r1")
  if (r2start<date & r2end> date) return("r2")
  if (r3start<date & r3end> date) return("r3")
  if (r4start<date & r4end> date) return("r4")
  if (r5start<date & r5end> date) return("r5")
  if (r6start<date & r6end> date) return("r6")
  else return (NA)
}

#da_class<-lapply(da,classround) #not working 
#foo to classify da and fill in a vector 
dclass<-vector()
for (i in 1:length(da)){
  print(i)
  dclass[i]<-classround(da[i])
}
daclass_f<-as.factor(dclass) #curse class as.factor 
summary(daclass_f) #summary 

#summary(daclass_f[newdata$idcol=="E049"]) # a test 
for( i in 1:length(unique(newdata$idcol))){
  id<-unique(newdata$idcol)[i]
  print(id)
  print(summary(daclass_f[newdata$idcol==id]))
}

#######################################################################################################################
######################################################################################################################
#select an ID of GPS collar that at least recorded sufficient information for an animal group in every season.
#create a new variable called Rounds in the data frame to identify relocation points based on the UTC_Date  

newdata$Rounds <- daclass_f
head(newdata)
table(newdata$idcol, newdata$Rounds)

#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E034' & newdata$Rounds=='r1'),] #24246 obs.
newdata.r2 <- newdata[which(newdata$idcol=='E037' & newdata$Rounds=='r2'),]#44613
newdata.r3 <- newdata[which(newdata$idcol=='e037' & newdata$Rounds=='r3'),]# 39966
newdata.r4 <- newdata[which(newdata$idcol=='e037' & newdata$Rounds=='r4'),]#47475
newdata.r5 <- newdata[which(newdata$idcol=='e037' & newdata$Rounds=='r5'),]#39186
newdata.r6 <- newdata[which(newdata$idcol=='e037' & newdata$Rounds=='r6'),]# 52728


hhid.184.herd.f <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.184.herd.f)
plot(hhid.184.herd.f$Longitude, hhid.184.herd.f$Latitude, col=hhid.184.herd.f$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  
hhid.184.herd.f.sub<- subset(hhid.184.herd.f, hhid.184.herd.f$Longitude>38.5 & hhid.184.herd.f$Longitude<39.5 & hhid.184.herd.f$Latitude>4.0 & hhid.184.herd.f$Latitude<4.5)
plot(hhid.184.herd.f.sub$Longitude, hhid.184.herd.f.sub$Latitude, col=hhid.184.herd.f.sub$Rounds)

#create a uniquely new ID for these collars.  
hhid.184.herd.f.sub$idcol_new <- as.factor("c006")

#ensure equal length for key variables 
head(hhid.184.herd.f.sub)
length(hhid.184.herd.f.sub$idcol_new) 
length(hhid.184.herd.f.sub$Longitude)
length(hhid.184.herd.f.sub$Latitude) 
length(hhid.184.herd.f.sub$UTC_Date) 
length(hhid.184.herd.f.sub$UTC_Time) 



saveRDS(hhid.184.herd.f.sub, file = "c006.hhid.184.fora.subherd.rds")

newdata <- readRDS("c006.hhid.184.fora.subherd.rds")
plot(newdata$Longitude, newdata$Latitude, col=newdata$Rounds)
summary(newdata$idcol_new)


#################################################################################################################
#################################################################################################################

###hhid 206, name Kanee Karayuu Didaa   

################################################################################################################
################################################################################################################
#import csv files 

list.files()
list_files <- list.files()
length <- length(list_files)

##write cattleTracks to R data frame data type
###identiy a group and animal
start <- proc.time()
newdata<-data.frame()# pre-allocate some variable space for our output data
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

head(newdata)
#p.730 <- newdata
saveRDS(newdata,file="p.206.allrounds.data.raw.rds") 

sum(is.na(newdata$Longitude)) #a total of 0 nas
sum(is.na(newdata$Latitude)) #a total of 0 nas 
sum(is.na(newdata$Elevation)) #a total of 0 nas
sum(is.na(newdata$UTC_Date))
sum(is.na(newdata$UTC_Time))

#choose to run these lines only to select complete data without nas
# n.data <- readRDS("p.155.allrounds.data.raw.rds") 
# n.data.2 <- n.data[!is.na(n.data$Longitude),]
# head(n.data.2)
# sum(is.na(n.data.2$Longitude)) #a total of 0 nas
# sum(is.na(n.data.2$Latitude)) #a total of 0 nas 
# sum(is.na(n.data.2$Elevation)) 
# sum(is.na(n.data.2$UTC_Date))
# sum(is.na(n.data.2$UTC_Time))

#new_DF <- DF[is.na(DF$Var),] #an example



plot(newdata$Longitude, newdata$Latitude, col=newdata$idcol) #pastoralist had both temporary and permanent place since the inception of study.

da<-strptime(paste(newdata$UTC_Date,newdata$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #mario 

#r1
r1start<-as.POSIXct("2011-08-01 00:00:00 CEST")
r1end<-as.POSIXct("2012-03-31 00:00:00 CEST")
#r2
r2start<-as.POSIXct("2012-09-01 00:00:00 CEST")
r2end<-as.POSIXct("2013-02-27 00:00:00 CEST")
#r3
r3start<-as.POSIXct("2013-09-01 00:00:00 CEST")
r3end<-as.POSIXct("2014-02-14 00:00:00 CEST")
#r4
r4start<-as.POSIXct("2014-02-15 00:00:00 CEST")
r4end<-as.POSIXct("2014-08-31 00:00:00 CEST")
#r5
r5start<-as.POSIXct("2014-09-01 00:00:00 CEST")
r5end<-as.POSIXct("2015-02-12 00:00:00 CEST")
#r6
r6start<-as.POSIXct("2015-02-13 00:00:00 CEST")#not 15th feb but 13th feb, and contains about 292 NAs
r6end<-as.POSIXct("2015-09-30 00:00:00 CEST") #changed to sep 30 from aug 31

da[2]<da[1]
classround<-function(date){
  if (is.na(date)) return (NA)
  if (r1start<date & r1end> date) return("r1")
  if (r2start<date & r2end> date) return("r2")
  if (r3start<date & r3end> date) return("r3")
  if (r4start<date & r4end> date) return("r4")
  if (r5start<date & r5end> date) return("r5")
  if (r6start<date & r6end> date) return("r6")
  else return (NA)
}

#da_class<-lapply(da,classround) #not working 
#foo to classify da and fill in a vector 
dclass<-vector()
for (i in 1:length(da)){
  print(i)
  dclass[i]<-classround(da[i])
}
daclass_f<-as.factor(dclass) #curse class as.factor 
summary(daclass_f) #summary 

#summary(daclass_f[newdata$idcol=="E049"]) # a test 
for( i in 1:length(unique(newdata$idcol))){
  id<-unique(newdata$idcol)[i]
  print(id)
  print(summary(daclass_f[newdata$idcol==id]))
}

#######################################################################################################################
######################################################################################################################
#select an ID of GPS collar that at least recorded sufficient information for an animal group in every season.
#create a new variable called Rounds in the data frame to identify relocation points based on the UTC_Date  

newdata$Rounds <- daclass_f
head(newdata)
table(newdata$idcol, newdata$Rounds)

#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E028' & newdata$Rounds=='r1'),] #29789 obs.
newdata.r2 <- newdata[which(newdata$idcol=='E007' & newdata$Rounds=='r2'),]#44596
newdata.r3 <- newdata[which(newdata$idcol=='e007' & newdata$Rounds=='r3'),]# 38228
newdata.r4 <- newdata[which(newdata$idcol=='e007' & newdata$Rounds=='r4'),]#47141
newdata.r5 <- newdata[which(newdata$idcol=='e007' & newdata$Rounds=='r5'),]#39365
newdata.r6 <- newdata[which(newdata$idcol=='r067' & newdata$Rounds=='r6'),]# 31164


hhid.206.herd.f <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.206.herd.f)
plot(hhid.206.herd.f$Longitude, hhid.206.herd.f$Latitude, col=hhid.206.herd.f$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  
hhid.206.herd.f.sub<- subset(hhid.206.herd.f, hhid.206.herd.f$Longitude>38.45 & hhid.206.herd.f$Longitude<38.9 & hhid.206.herd.f$Latitude>4.0 & hhid.206.herd.f$Latitude<5.0)
plot(hhid.206.herd.f.sub$Longitude, hhid.206.herd.f.sub$Latitude, col=hhid.206.herd.f.sub$Rounds)

#create a uniquely new ID for these collars.  
hhid.206.herd.f.sub$idcol_new <- as.factor("c007")

#ensure equal length for key variables 
head(hhid.206.herd.f.sub)
length(hhid.206.herd.f.sub$idcol_new) 
length(hhid.206.herd.f.sub$Longitude)
length(hhid.206.herd.f.sub$Latitude) 
length(hhid.206.herd.f.sub$UTC_Date) 
length(hhid.206.herd.f.sub$UTC_Time) 



saveRDS(hhid.206.herd.f.sub, file = "c007.hhid.206.fora.subherd.rds")

newdata <- readRDS("c007.hhid.206.fora.subherd.rds")
plot(newdata$Longitude, newdata$Latitude, col=newdata$Rounds)
summary(newdata$idcol_new)

#################################################################################################################
#################################################################################################################

###hhid 468, name Jilgaa Raqo Areero  

################################################################################################################
################################################################################################################
#import csv files 

list.files()
list_files <- list.files()
length <- length(list_files)

##write cattleTracks to R data frame data type
###identiy a group and animal
start <- proc.time()
newdata<-data.frame()# pre-allocate some variable space for our output data
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

head(newdata)
#p.730 <- newdata
saveRDS(newdata,file="p.468.allrounds.data.raw.rds") 
newdata <- readRDS("p.468.allrounds.data.raw.rds")

sum(is.na(newdata$Longitude)) #a total of 0 nas
sum(is.na(newdata$Latitude)) #a total of 0 nas 
sum(is.na(newdata$Elevation)) #a total of 0 nas
sum(is.na(newdata$UTC_Date))
sum(is.na(newdata$UTC_Time))

#choose to run these lines only to select complete data without nas
# n.data <- readRDS("p.155.allrounds.data.raw.rds") 
# n.data.2 <- n.data[!is.na(n.data$Longitude),]
# head(n.data.2)
# sum(is.na(n.data.2$Longitude)) #a total of 0 nas
# sum(is.na(n.data.2$Latitude)) #a total of 0 nas 
# sum(is.na(n.data.2$Elevation)) 
# sum(is.na(n.data.2$UTC_Date))
# sum(is.na(n.data.2$UTC_Time))

#new_DF <- DF[is.na(DF$Var),] #an example



plot(newdata$Longitude, newdata$Latitude, col=newdata$idcol) #pastoralist had both temporary and permanent place since the inception of study.

da<-strptime(paste(newdata$UTC_Date,newdata$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #mario 

#r1
r1start<-as.POSIXct("2011-08-01 00:00:00 CEST")
r1end<-as.POSIXct("2012-03-31 00:00:00 CEST")
#r2
r2start<-as.POSIXct("2012-09-01 00:00:00 CEST")
r2end<-as.POSIXct("2013-02-27 00:00:00 CEST")
#r3
r3start<-as.POSIXct("2013-09-01 00:00:00 CEST")
r3end<-as.POSIXct("2014-02-14 00:00:00 CEST")
#r4
r4start<-as.POSIXct("2014-02-15 00:00:00 CEST")
r4end<-as.POSIXct("2014-08-31 00:00:00 CEST")
#r5
r5start<-as.POSIXct("2014-09-01 00:00:00 CEST")
r5end<-as.POSIXct("2015-02-12 00:00:00 CEST")
#r6
r6start<-as.POSIXct("2015-02-13 00:00:00 CEST")#not 15th feb but 13th feb, and contains about 292 NAs
r6end<-as.POSIXct("2015-09-30 00:00:00 CEST") #changed to sep 30 from aug 31

da[2]<da[1]
classround<-function(date){
  if (is.na(date)) return (NA)
  if (r1start<date & r1end> date) return("r1")
  if (r2start<date & r2end> date) return("r2")
  if (r3start<date & r3end> date) return("r3")
  if (r4start<date & r4end> date) return("r4")
  if (r5start<date & r5end> date) return("r5")
  if (r6start<date & r6end> date) return("r6")
  else return (NA)
}

#da_class<-lapply(da,classround) #not working 
#foo to classify da and fill in a vector 
dclass<-vector()
for (i in 1:length(da)){
  print(i)
  dclass[i]<-classround(da[i])
}
daclass_f<-as.factor(dclass) #curse class as.factor 
summary(daclass_f) #summary 

#summary(daclass_f[newdata$idcol=="E049"]) # a test 
for( i in 1:length(unique(newdata$idcol))){
  id<-unique(newdata$idcol)[i]
  print(id)
  print(summary(daclass_f[newdata$idcol==id]))
}

#######################################################################################################################
######################################################################################################################
#select an ID of GPS collar that at least recorded sufficient information for an animal group in every season.
#create a new variable called Rounds in the data frame to identify relocation points based on the UTC_Date  

newdata$Rounds <- daclass_f
head(newdata)
table(newdata$idcol, newdata$Rounds) #r5, E107 belongs to a different hhid.. confirm whether this claim is justified  

#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E019' & newdata$Rounds=='r1'),] #59667 obs.
newdata.r2 <- newdata[which(newdata$idcol=='E020' & newdata$Rounds=='r2'),]#46350
#newdata.r3 <- newdata[which(newdata$idcol=='e020' & newdata$Rounds=='r3'),]# 4424 #probably, a collar was deployed to household in dillo 
newdata.r4 <- newdata[which(newdata$idcol=='e107' & newdata$Rounds=='r4'),]#36322
newdata.r5 <- newdata[which(newdata$idcol=='e110' & newdata$Rounds=='r5'),]#42668
newdata.r6 <- newdata[which(newdata$idcol=='e110' & newdata$Rounds=='r6'),]# 32153


hhid.468.herd.f <- rbind(newdata.r1, newdata.r2, newdata.r4, newdata.r5, newdata.r6)
head(hhid.468.herd.f)
plot(hhid.468.herd.f$Longitude, hhid.468.herd.f$Latitude, col=hhid.468.herd.f$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning. all households in sarite subset at lon 37.6, especially round1 data  
hhid.468.herd.f.sub<- subset(hhid.468.herd.f, hhid.468.herd.f$Longitude>37.52 & hhid.468.herd.f$Longitude<37.72 & hhid.468.herd.f$Latitude>4.0 & hhid.468.herd.f$Latitude<5.0)
plot(hhid.468.herd.f.sub$Longitude, hhid.468.herd.f.sub$Latitude, col=hhid.468.herd.f.sub$Rounds)

#create a uniquely new ID for these collars.  
hhid.468.herd.f.sub$idcol_new <- as.factor("c008")

#ensure equal length for key variables 
head(hhid.468.herd.f.sub)
length(hhid.468.herd.f.sub$idcol_new) 
length(hhid.468.herd.f.sub$Longitude)
length(hhid.468.herd.f.sub$Latitude) 
length(hhid.468.herd.f.sub$UTC_Date) 
length(hhid.468.herd.f.sub$UTC_Time) 



saveRDS(hhid.468.herd.f.sub, file = "c008.hhid.468.wora.subherd.rds")
newdata <- readRDS("c008.hhid.468.wora.subherd.rds")


#################################################################################################################
#################################################################################################################

###hhid 540, name Jarso Qancora Mamo 

################################################################################################################
################################################################################################################
#import csv files 

list.files()
list_files <- list.files()
length <- length(list_files)

##write cattleTracks to R data frame data type
###identiy a group and animal
start <- proc.time()
newdata<-data.frame()# pre-allocate some variable space for our output data
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

head(newdata)
#p.730 <- newdata
saveRDS(newdata,file="p.540.allrounds.data.raw.rds") 
#newdata <-readRDS("p.540.allrounds.data.raw.rds") 

sum(is.na(newdata$Longitude)) #a total of 0 nas
sum(is.na(newdata$Latitude)) #a total of 0 nas 
sum(is.na(newdata$Elevation)) #a total of 0 nas
sum(is.na(newdata$UTC_Date))
sum(is.na(newdata$UTC_Time))

# #choose to run these lines only to select complete data without nas
# n.data <- readRDS("p.540.allrounds.data.raw.rds")
# n.data.2 <- n.data[!is.na(n.data$Longitude),]
# head(n.data.2)
# sum(is.na(n.data.2$Longitude)) #a total of 0 nas
# sum(is.na(n.data.2$Latitude)) #a total of 0 nas
# sum(is.na(n.data.2$Elevation))
# sum(is.na(n.data.2$UTC_Date))
# sum(is.na(n.data.2$UTC_Time))

#new_DF <- DF[is.na(DF$Var),] #an example
#newdata <- n.data.2


plot(newdata$Longitude, newdata$Latitude, col=newdata$idcol) #pastoralist had both temporary and permanent place since the inception of study.

da<-strptime(paste(newdata$UTC_Date,newdata$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #mario 

#r1
r1start<-as.POSIXct("2011-08-01 00:00:00 CEST")
r1end<-as.POSIXct("2012-03-31 00:00:00 CEST")
#r2
r2start<-as.POSIXct("2012-09-01 00:00:00 CEST")
r2end<-as.POSIXct("2013-02-27 00:00:00 CEST")
#r3
r3start<-as.POSIXct("2013-09-01 00:00:00 CEST")
r3end<-as.POSIXct("2014-02-14 00:00:00 CEST")
#r4
r4start<-as.POSIXct("2014-02-15 00:00:00 CEST")
r4end<-as.POSIXct("2014-08-31 00:00:00 CEST")
#r5
r5start<-as.POSIXct("2014-09-01 00:00:00 CEST")
r5end<-as.POSIXct("2015-02-12 00:00:00 CEST")
#r6
r6start<-as.POSIXct("2015-02-13 00:00:00 CEST")#not 15th feb but 13th feb, and contains about 292 NAs
r6end<-as.POSIXct("2015-09-30 00:00:00 CEST") #changed to sep 30 from aug 31

da[2]<da[1]
classround<-function(date){
  if (is.na(date)) return (NA)
  if (r1start<date & r1end> date) return("r1")
  if (r2start<date & r2end> date) return("r2")
  if (r3start<date & r3end> date) return("r3")
  if (r4start<date & r4end> date) return("r4")
  if (r5start<date & r5end> date) return("r5")
  if (r6start<date & r6end> date) return("r6")
  else return (NA)
}

#da_class<-lapply(da,classround) #not working 
#foo to classify da and fill in a vector 
dclass<-vector()
for (i in 1:length(da)){
  print(i)
  dclass[i]<-classround(da[i])
}
daclass_f<-as.factor(dclass) #curse class as.factor 
summary(daclass_f) #summary 

#summary(daclass_f[newdata$idcol=="E049"]) # a test 
for( i in 1:length(unique(newdata$idcol))){
  id<-unique(newdata$idcol)[i]
  print(id)
  print(summary(daclass_f[newdata$idcol==id]))
}

#######################################################################################################################
######################################################################################################################
#select an ID of GPS collar that at least recorded sufficient information for an animal group in every season.
#create a new variable called Rounds in the data frame to identify relocation points based on the UTC_Date  

newdata$Rounds <- daclass_f
head(newdata)
table(newdata$idcol, newdata$Rounds) #  

#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E022' & newdata$Rounds=='r1'),] #54278 obs.
newdata.r2 <- newdata[which(newdata$idcol=='E014' & newdata$Rounds=='r2'),]#46441
newdata.r3 <- newdata[which(newdata$idcol=='e014' & newdata$Rounds=='r3'),]# 40409
newdata.r4 <- newdata[which(newdata$idcol=='e024' & newdata$Rounds=='r4'),]#46723
newdata.r5 <- newdata[which(newdata$idcol=='e024' & newdata$Rounds=='r5'),]#7120
newdata.r6 <- newdata[which(newdata$idcol=='e014' & newdata$Rounds=='r6'),]# 50977


hhid.540.herd.f <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.540.herd.f)
plot(hhid.540.herd.f$Longitude, hhid.540.herd.f$Latitude, col=hhid.540.herd.f$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  
hhid.540.herd.f.sub<- subset(hhid.540.herd.f, hhid.540.herd.f$Longitude>37.4 & hhid.540.herd.f$Longitude<37.75 & hhid.540.herd.f$Latitude>4.5 & hhid.540.herd.f$Latitude<5.1)
plot(hhid.540.herd.f.sub$Longitude, hhid.540.herd.f.sub$Latitude, col=hhid.540.herd.f.sub$Rounds)

#create a uniquely new ID for these collars.  
hhid.540.herd.f.sub$idcol_new <- as.factor("c009")

#ensure equal length for key variables 
head(hhid.540.herd.f.sub)
length(hhid.540.herd.f.sub$idcol_new) 
length(hhid.540.herd.f.sub$Longitude)
length(hhid.540.herd.f.sub$Latitude) 
length(hhid.540.herd.f.sub$UTC_Date) 
length(hhid.540.herd.f.sub$UTC_Time) 



saveRDS(hhid.540.herd.f.sub, file = "c009.hhid.540.wora.subherd.rds")


#################################################################################################################
#################################################################################################################

###hhid 566, name Elema Galma 

################################################################################################################
################################################################################################################
#import csv files 

list.files()
list_files <- list.files()
length <- length(list_files)

##write cattleTracks to R data frame data type
###identiy a group and animal
start <- proc.time()
newdata<-data.frame()# pre-allocate some variable space for our output data
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

head(newdata)
#p.730 <- newdata
saveRDS(newdata,file="p.566.allrounds.data.raw.rds")
newdata <- readRDS("p.566.allrounds.data.raw.rds")

sum(is.na(newdata$Longitude)) #a total of 0 nas
sum(is.na(newdata$Latitude)) #a total of 0 nas 
sum(is.na(newdata$Elevation)) #a total of 0 nas
sum(is.na(newdata$UTC_Date))
sum(is.na(newdata$UTC_Time))

# #choose to run these lines only to select complete data without nas
# n.data <- readRDS("p.540.allrounds.data.raw.rds")
# n.data.2 <- n.data[!is.na(n.data$Longitude),]
# head(n.data.2)
# sum(is.na(n.data.2$Longitude)) #a total of 0 nas
# sum(is.na(n.data.2$Latitude)) #a total of 0 nas
# sum(is.na(n.data.2$Elevation))
# sum(is.na(n.data.2$UTC_Date))
# sum(is.na(n.data.2$UTC_Time))
# 
# #new_DF <- DF[is.na(DF$Var),] #an example
# newdata <- n.data.2


plot(newdata$Longitude, newdata$Latitude, col=newdata$idcol) #pastoralist had both temporary and permanent place since the inception of study.

da<-strptime(paste(newdata$UTC_Date,newdata$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #mario 

#r1
r1start<-as.POSIXct("2011-08-01 00:00:00 CEST")
r1end<-as.POSIXct("2012-03-31 00:00:00 CEST")
#r2
r2start<-as.POSIXct("2012-09-01 00:00:00 CEST")
r2end<-as.POSIXct("2013-02-27 00:00:00 CEST")
#r3
r3start<-as.POSIXct("2013-09-01 00:00:00 CEST")
r3end<-as.POSIXct("2014-02-14 00:00:00 CEST")
#r4
r4start<-as.POSIXct("2014-02-15 00:00:00 CEST")
r4end<-as.POSIXct("2014-08-31 00:00:00 CEST")
#r5
r5start<-as.POSIXct("2014-09-01 00:00:00 CEST")
r5end<-as.POSIXct("2015-02-12 00:00:00 CEST")
#r6
r6start<-as.POSIXct("2015-02-13 00:00:00 CEST")#not 15th feb but 13th feb, and contains about 292 NAs
r6end<-as.POSIXct("2015-09-30 00:00:00 CEST") #changed to sep 30 from aug 31

da[2]<da[1]
classround<-function(date){
  if (is.na(date)) return (NA)
  if (r1start<date & r1end> date) return("r1")
  if (r2start<date & r2end> date) return("r2")
  if (r3start<date & r3end> date) return("r3")
  if (r4start<date & r4end> date) return("r4")
  if (r5start<date & r5end> date) return("r5")
  if (r6start<date & r6end> date) return("r6")
  else return (NA)
}

#da_class<-lapply(da,classround) #not working 
#foo to classify da and fill in a vector 
dclass<-vector()
for (i in 1:length(da)){
  print(i)
  dclass[i]<-classround(da[i])
}
daclass_f<-as.factor(dclass) #curse class as.factor 
summary(daclass_f) #summary 

#summary(daclass_f[newdata$idcol=="E049"]) # a test 
for( i in 1:length(unique(newdata$idcol))){
  id<-unique(newdata$idcol)[i]
  print(id)
  print(summary(daclass_f[newdata$idcol==id]))
}

#######################################################################################################################
######################################################################################################################
#select an ID of GPS collar that at least recorded sufficient information for an animal group in every season.
#create a new variable called Rounds in the data frame to identify relocation points based on the UTC_Date  

newdata$Rounds <- daclass_f
head(newdata)
table(newdata$idcol, newdata$Rounds) #  

#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E018' & newdata$Rounds=='r1'),] #58759 obs.
newdata.r2 <- newdata[which(newdata$idcol=='E016' & newdata$Rounds=='r2'),]#
newdata.r3 <- newdata[which(newdata$idcol=='e061' & newdata$Rounds=='r3'),]# 
newdata.r4 <- newdata[which(newdata$idcol=='e101' & newdata$Rounds=='r4'),]#
newdata.r5 <- newdata[which(newdata$idcol=='r064' & newdata$Rounds=='r5'),]#
newdata.r6 <- newdata[which(newdata$idcol=='e047' & newdata$Rounds=='r6'),]# 


hhid.566.herd.w <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.566.herd.w)
plot(hhid.566.herd.w$Longitude, hhid.566.herd.w$Latitude, col=hhid.566.herd.w$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  
hhid.566.herd.w.sub<- subset(hhid.566.herd.w, hhid.566.herd.w$Longitude>37.53 & hhid.566.herd.w$Longitude<37.8 & hhid.566.herd.w$Latitude>4.8 & hhid.566.herd.w$Latitude<5.0)
plot(hhid.566.herd.w.sub$Longitude, hhid.566.herd.w.sub$Latitude, col=hhid.566.herd.w.sub$Rounds)

#create a uniquely new ID for these collars.  
hhid.566.herd.w.sub$idcol_new <- as.factor("c010")
summary(hhid.566.herd.w.sub$idcol_new)

#ensure equal length for key variables 
head(hhid.566.herd.w.sub)
length(hhid.566.herd.w.sub$idcol_new) 
length(hhid.566.herd.w.sub$Longitude)
length(hhid.566.herd.w.sub$Latitude) 
length(hhid.566.herd.w.sub$UTC_Date) 
length(hhid.566.herd.w.sub$UTC_Time) 



saveRDS(hhid.566.herd.w.sub, file = "c010.w.hhid.566.wora.subherd.rds") #wora sub herd


#compute a trajectory for a sub herd fora belonging to household hhid 566 
#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E017' & newdata$Rounds=='r1'),] #58759 obs.
newdata.r2 <- newdata[which(newdata$idcol=='E018' & newdata$Rounds=='r2'),]#
newdata.r3 <- newdata[which(newdata$idcol=='e062' & newdata$Rounds=='r3'),]# 
newdata.r4 <- newdata[which(newdata$idcol=='r064' & newdata$Rounds=='r4'),]#
newdata.r5 <- newdata[which(newdata$idcol=='e101' & newdata$Rounds=='r5'),]#
newdata.r6 <- newdata[which(newdata$idcol=='e012' & newdata$Rounds=='r6'),]# 


hhid.566.herd.f <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.566.herd.f)
plot(hhid.566.herd.f$Longitude, hhid.566.herd.f$Latitude, col=hhid.566.herd.f$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  
hhid.566.herd.f.sub<- subset(hhid.566.herd.f, hhid.566.herd.f$Longitude>37.0 & hhid.566.herd.f$Longitude<37.72 & hhid.566.herd.f$Latitude>4.8 & hhid.566.herd.f$Latitude<5.0)
plot(hhid.566.herd.f.sub$Longitude, hhid.566.herd.f.sub$Latitude, col=hhid.566.herd.f.sub$Rounds)

#create a uniquely new ID for these collars.  
hhid.566.herd.f.sub$idcol_new <- as.factor("c010.f") #fora sub herd 
summary(hhid.566.herd.f.sub$idcol_new)

#ensure equal length for key variables 
head(hhid.566.herd.f.sub)
length(hhid.566.herd.f.sub$idcol_new) 
length(hhid.566.herd.f.sub$Longitude)
length(hhid.566.herd.f.sub$Latitude) 
length(hhid.566.herd.f.sub$UTC_Date) 
length(hhid.566.herd.f.sub$UTC_Time) 



saveRDS(hhid.566.herd.f.sub, file = "c010.hhid.566.wora.subherd.rds") #a fora sub herd 

#################################################################################################################
#################################################################################################################

###hhid 616, name Guyo Shiriti 

################################################################################################################
################################################################################################################
#import csv files 

list.files()
list_files <- list.files()
length <- length(list_files)

##write cattleTracks to R data frame data type
###identiy a group and animal
start <- proc.time()
newdata<-data.frame()# pre-allocate some variable space for our output data
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

head(newdata)
#p.730 <- newdata
saveRDS(newdata,file="p.616.allrounds.data.raw.rds")
newdata <- readRDS("p.616.allrounds.data.raw.rds")

sum(is.na(newdata$Longitude)) #a total of 0 nas
sum(is.na(newdata$Latitude)) #a total of 0 nas 
sum(is.na(newdata$Elevation)) #a total of 0 nas
sum(is.na(newdata$UTC_Date))
sum(is.na(newdata$UTC_Time))

# #choose to run these lines only to select complete data without nas
# n.data <- readRDS("p.540.allrounds.data.raw.rds")
# n.data.2 <- n.data[!is.na(n.data$Longitude),]
# head(n.data.2)
# sum(is.na(n.data.2$Longitude)) #a total of 0 nas
# sum(is.na(n.data.2$Latitude)) #a total of 0 nas
# sum(is.na(n.data.2$Elevation))
# sum(is.na(n.data.2$UTC_Date))
# sum(is.na(n.data.2$UTC_Time))
# 
# #new_DF <- DF[is.na(DF$Var),] #an example
# newdata <- n.data.2


plot(newdata$Longitude, newdata$Latitude, col=newdata$idcol) #pastoralist had both temporary and permanent place since the inception of study.

da<-strptime(paste(newdata$UTC_Date,newdata$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #mario 

#r1
r1start<-as.POSIXct("2011-08-01 00:00:00 CEST")
r1end<-as.POSIXct("2012-03-31 00:00:00 CEST")
#r2
r2start<-as.POSIXct("2012-09-01 00:00:00 CEST")
r2end<-as.POSIXct("2013-02-27 00:00:00 CEST")
#r3
r3start<-as.POSIXct("2013-09-01 00:00:00 CEST")
r3end<-as.POSIXct("2014-02-14 00:00:00 CEST")
#r4
r4start<-as.POSIXct("2014-02-15 00:00:00 CEST")
r4end<-as.POSIXct("2014-08-31 00:00:00 CEST")
#r5
r5start<-as.POSIXct("2014-09-01 00:00:00 CEST")
r5end<-as.POSIXct("2015-02-12 00:00:00 CEST")
#r6
r6start<-as.POSIXct("2015-02-13 00:00:00 CEST")#not 15th feb but 13th feb, and contains about 292 NAs
r6end<-as.POSIXct("2015-09-30 00:00:00 CEST") #changed to sep 30 from aug 31

da[2]<da[1]
classround<-function(date){
  if (is.na(date)) return (NA)
  if (r1start<date & r1end> date) return("r1")
  if (r2start<date & r2end> date) return("r2")
  if (r3start<date & r3end> date) return("r3")
  if (r4start<date & r4end> date) return("r4")
  if (r5start<date & r5end> date) return("r5")
  if (r6start<date & r6end> date) return("r6")
  else return (NA)
}

#da_class<-lapply(da,classround) #not working 
#foo to classify da and fill in a vector 
dclass<-vector()
for (i in 1:length(da)){
  print(i)
  dclass[i]<-classround(da[i])
}
daclass_f<-as.factor(dclass) #curse class as.factor 
summary(daclass_f) #summary 

#summary(daclass_f[newdata$idcol=="E049"]) # a test 
for( i in 1:length(unique(newdata$idcol))){
  id<-unique(newdata$idcol)[i]
  print(id)
  print(summary(daclass_f[newdata$idcol==id]))
}

#######################################################################################################################
######################################################################################################################
#select an ID of GPS collar that at least recorded sufficient information for an animal group in every season.
#create a new variable called Rounds in the data frame to identify relocation points based on the UTC_Date  

newdata$Rounds <- daclass_f
head(newdata)
table(newdata$idcol, newdata$Rounds) #r5, E107 belongs to a different hhid.. confirm whether this claim is justified  

#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E013' & newdata$Rounds=='r1'),] #obs.
newdata.r2 <- newdata[which(newdata$idcol=='E057' & newdata$Rounds=='r2'),]#
newdata.r3 <- newdata[which(newdata$idcol=='e057' & newdata$Rounds=='r3'),]# 
newdata.r4 <- newdata[which(newdata$idcol=='e057' & newdata$Rounds=='r4'),]#
newdata.r5 <- newdata[which(newdata$idcol=='e063' & newdata$Rounds=='r5'),]#
newdata.r6 <- newdata[which(newdata$idcol=='e063' & newdata$Rounds=='r6'),]# 

##special cleaning of round 1
 
# newdata.r1.1 <- subset(newdata.r1, newdata.r1$Longitude>37.6 & newdata.r1$Longitude<37.75 & newdata.r1$Latitude>3.5 & newdata.r1$Latitude<5.0)
# plot(newdata.r1.1$Longitude, newdata.r1.1$Latitude, col=newdata.r1.1$Rounds)
# newdata.r1 <- newdata.r1.1

hhid.616.herd.f <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.616.herd.f)
plot(hhid.616.herd.f$Longitude, hhid.616.herd.f$Latitude, col=hhid.616.herd.f$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  take lon greater than 37.6 for all households in sarite except for guyo shiriti
hhid.616.herd.f.sub<- subset(hhid.616.herd.f, hhid.616.herd.f$Longitude>37.2 & hhid.616.herd.f$Longitude<37.9 & hhid.616.herd.f$Latitude>3.5 & hhid.616.herd.f$Latitude<5.4)
plot(hhid.616.herd.f.sub$Longitude, hhid.616.herd.f.sub$Latitude, col=hhid.616.herd.f.sub$Rounds)



#create a uniquely new ID for these collars.  
hhid.616.herd.f.sub$idcol_new <- as.factor("c011.w") #add a suffix to as.factor in this line to identify a multiple herds 

#ensure equal length for key variables 
head(hhid.616.herd.f.sub)
length(hhid.616.herd.f.sub$idcol_new) 
length(hhid.616.herd.f.sub$Longitude)
length(hhid.616.herd.f.sub$Latitude) 
length(hhid.616.herd.f.sub$UTC_Date) 
length(hhid.616.herd.f.sub$UTC_Time) 



saveRDS(hhid.616.herd.f.sub, file = "c011.hhid.616.wora.subherd.rds") #wara sub herd 

###create a dataframe for a fora sub herd 
#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E014' & newdata$Rounds=='r1'),] #obs.
newdata.r2 <- newdata[which(newdata$idcol=='E056' & newdata$Rounds=='r2'),]#
newdata.r3 <- newdata[which(newdata$idcol=='e063' & newdata$Rounds=='r3'),]# 
newdata.r4 <- newdata[which(newdata$idcol=='e063' & newdata$Rounds=='r4'),]#
newdata.r5 <- newdata[which(newdata$idcol=='e106' & newdata$Rounds=='r5'),]#
newdata.r6 <- newdata[which(newdata$idcol=='e106' & newdata$Rounds=='r6'),]# 

##special cleaning of round 1

# newdata.r1.1 <- subset(newdata.r1, newdata.r1$Longitude>37.6 & newdata.r1$Longitude<37.75 & newdata.r1$Latitude>3.5 & newdata.r1$Latitude<5.0)
# plot(newdata.r1.1$Longitude, newdata.r1.1$Latitude, col=newdata.r1.1$Rounds)
# newdata.r1 <- newdata.r1.1

hhid.616.herd.w <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.616.herd.w)
plot(hhid.616.herd.w$Longitude, hhid.616.herd.w$Latitude, col=hhid.616.herd.w$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  take lon greater than 37.6 for all households in sarite except for guyo shiriti
hhid.616.herd.w.sub<- subset(hhid.616.herd.w, hhid.616.herd.w$Longitude>37.5 & hhid.616.herd.w$Longitude<37.75 & hhid.616.herd.w$Latitude>3.5 & hhid.616.herd.w$Latitude<5.0)
plot(hhid.616.herd.w.sub$Longitude, hhid.616.herd.w.sub$Latitude, col=hhid.616.herd.w.sub$Rounds)



#create a uniquely new ID for these collars.  
hhid.616.herd.w.sub$idcol_new <- as.factor("c011") #add a suffix to as.factor in this line to identify a multiple herds 

#ensure equal length for key variables 
head(hhid.616.herd.w.sub)
length(hhid.616.herd.w.sub$idcol_new) 
length(hhid.616.herd.w.sub$Longitude)
length(hhid.616.herd.w.sub$Latitude) 
length(hhid.616.herd.w.sub$UTC_Date) 
length(hhid.616.herd.w.sub$UTC_Time) 



saveRDS(hhid.616.herd.w.sub, file = "c011.w.hhid.616.wora.subherd.rds") #wara sub herd 



########arero

#################################################################################################################
#################################################################################################################

###hhid 34, name Tache Kote Abbakiyo

################################################################################################################
################################################################################################################
#import csv files 

list.files()
list_files <- list.files()
length <- length(list_files)

##write cattleTracks to R data frame data type
###identiy a group and animal
start <- proc.time()
newdata<-data.frame()# pre-allocate some variable space for our output data
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

head(newdata)
#p.730 <- newdata
saveRDS(newdata,file="p.034.allrounds.data.raw.rds") 

sum(is.na(newdata$Longitude)) #a total of 0 nas
sum(is.na(newdata$Latitude)) #a total of 0 nas 
sum(is.na(newdata$Elevation)) #a total of 0 nas
sum(is.na(newdata$UTC_Date))
sum(is.na(newdata$UTC_Time))

# #choose to run these lines only to select complete data without nas
# n.data <- readRDS("p.540.allrounds.data.raw.rds")
# n.data.2 <- n.data[!is.na(n.data$Longitude),]
# head(n.data.2)
# sum(is.na(n.data.2$Longitude)) #a total of 0 nas
# sum(is.na(n.data.2$Latitude)) #a total of 0 nas
# sum(is.na(n.data.2$Elevation))
# sum(is.na(n.data.2$UTC_Date))
# sum(is.na(n.data.2$UTC_Time))
# 
# #new_DF <- DF[is.na(DF$Var),] #an example
# newdata <- n.data.2


plot(newdata$Longitude, newdata$Latitude, col=newdata$idcol) #pastoralist had both temporary and permanent place since the inception of study.

da<-strptime(paste(newdata$UTC_Date,newdata$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #mario 

#r1
r1start<-as.POSIXct("2011-08-01 00:00:00 CEST")
r1end<-as.POSIXct("2012-03-31 00:00:00 CEST")
#r2
r2start<-as.POSIXct("2012-09-01 00:00:00 CEST")
r2end<-as.POSIXct("2013-02-27 00:00:00 CEST")
#r3
r3start<-as.POSIXct("2013-09-01 00:00:00 CEST")
r3end<-as.POSIXct("2014-02-14 00:00:00 CEST")
#r4
r4start<-as.POSIXct("2014-02-15 00:00:00 CEST")
r4end<-as.POSIXct("2014-08-31 00:00:00 CEST")
#r5
r5start<-as.POSIXct("2014-09-01 00:00:00 CEST")
r5end<-as.POSIXct("2015-02-12 00:00:00 CEST")
#r6
r6start<-as.POSIXct("2015-02-13 00:00:00 CEST")#not 15th feb but 13th feb, and contains about 292 NAs
r6end<-as.POSIXct("2015-09-30 00:00:00 CEST") #changed to sep 30 from aug 31

da[2]<da[1]
classround<-function(date){
  if (is.na(date)) return (NA)
  if (r1start<date & r1end> date) return("r1")
  if (r2start<date & r2end> date) return("r2")
  if (r3start<date & r3end> date) return("r3")
  if (r4start<date & r4end> date) return("r4")
  if (r5start<date & r5end> date) return("r5")
  if (r6start<date & r6end> date) return("r6")
  else return (NA)
}

#da_class<-lapply(da,classround) #not working 
#foo to classify da and fill in a vector 
dclass<-vector()
for (i in 1:length(da)){
  print(i)
  dclass[i]<-classround(da[i])
}
daclass_f<-as.factor(dclass) #curse class as.factor 
summary(daclass_f) #summary 

#summary(daclass_f[newdata$idcol=="E049"]) # a test 
for( i in 1:length(unique(newdata$idcol))){
  id<-unique(newdata$idcol)[i]
  print(id)
  print(summary(daclass_f[newdata$idcol==id]))
}

#######################################################################################################################
######################################################################################################################
#select an ID of GPS collar that at least recorded sufficient information for an animal group in every season.
#create a new variable called Rounds in the data frame to identify relocation points based on the UTC_Date  

newdata$Rounds <- daclass_f
head(newdata)
table(newdata$idcol, newdata$Rounds) #r5, E107 belongs to a different hhid.. confirm whether this claim is justified  

#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E037' & newdata$Rounds=='r1'),] #obs.
newdata.r2 <- newdata[which(newdata$idcol=='E038' & newdata$Rounds=='r2'),]#
newdata.r3 <- newdata[which(newdata$idcol=='e043' & newdata$Rounds=='r3'),]# 
newdata.r4 <- newdata[which(newdata$idcol=='e043' & newdata$Rounds=='r4'),]#
newdata.r5 <- newdata[which(newdata$idcol=='e043' & newdata$Rounds=='r5'),]#
newdata.r6 <- newdata[which(newdata$idcol=='e038' & newdata$Rounds=='r6'),]# 

##special cleaning of round 1

# newdata.r1.1 <- subset(newdata.r1, newdata.r1$Longitude>37.6 & newdata.r1$Longitude<37.75 & newdata.r1$Latitude>3.5 & newdata.r1$Latitude<5.0)
# plot(newdata.r1.1$Longitude, newdata.r1.1$Latitude, col=newdata.r1.1$Rounds)
# newdata.r1 <- newdata.r1.1

hhid.034.herd.w <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.034.herd.w)
plot(hhid.034.herd.w$Longitude, hhid.034.herd.w$Latitude, col=hhid.034.herd.w$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  take lon greater than 37.6 for all households in sarite except for guyo shiriti
hhid.034.herd.f.sub<- subset(hhid.034.herd.f, hhid.034.herd.f$Longitude>38.0 & hhid.034.herd.f$Longitude<39.4 & hhid.034.herd.f$Latitude>4.0 & hhid.034.herd.f$Latitude<4.8)
plot(hhid.034.herd.f.sub$Longitude, hhid.034.herd.f.sub$Latitude, col=hhid.034.herd.f.sub$Rounds)

# cut at lon 39.1 for all the herd in arero

#create a uniquely new ID for these collars.  
hhid.034.herd.f.sub$idcol_new <- as.factor("c012") #add a suffix to as.factor in this line to identify a multiple herds 
#hhid.034.herd.w$idcol_new <- as.factor("c012.w")

#ensure equal length for key variables 
head(hhid.034.herd.f.sub)
length(hhid.034.herd.f.sub$idcol_new) 
length(hhid.034.herd.f.sub$Longitude)
length(hhid.034.herd.f.sub$Latitude) 
length(hhid.034.herd.f.sub$UTC_Date) 
length(hhid.034.herd.f.sub$UTC_Time) 



saveRDS(hhid.034.herd.f.sub, file = "c012.w.hhid.034.wora.subherd.rds")#? fora or wara herd, give a name to specify a herd. 

########arero


#################################################################################################################
#################################################################################################################

###hhid 46, name Wariyo Dido dhadacha 

################################################################################################################
################################################################################################################
#import csv files 

list.files()
list_files <- list.files()
length <- length(list_files)

##write cattleTracks to R data frame data type
###identiy a group and animal
start <- proc.time()
newdata<-data.frame()# pre-allocate some variable space for our output data
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

head(newdata)
#p.730 <- newdata
saveRDS(newdata,file="p.046.allrounds.data.raw.rds") #same trajectory data as 034, a suprise! 

sum(is.na(newdata$Longitude)) #a total of 0 nas
sum(is.na(newdata$Latitude)) #a total of 0 nas 
sum(is.na(newdata$Elevation)) #a total of 0 nas
sum(is.na(newdata$UTC_Date))
sum(is.na(newdata$UTC_Time))

# #choose to run these lines only to select complete data without nas
# n.data <- readRDS("p.540.allrounds.data.raw.rds")
# n.data.2 <- n.data[!is.na(n.data$Longitude),]
# head(n.data.2)
# sum(is.na(n.data.2$Longitude)) #a total of 0 nas
# sum(is.na(n.data.2$Latitude)) #a total of 0 nas
# sum(is.na(n.data.2$Elevation))
# sum(is.na(n.data.2$UTC_Date))
# sum(is.na(n.data.2$UTC_Time))
# 
# #new_DF <- DF[is.na(DF$Var),] #an example
# newdata <- n.data.2


plot(newdata$Longitude, newdata$Latitude, col=newdata$idcol) #pastoralist had both temporary and permanent place since the inception of study.

da<-strptime(paste(newdata$UTC_Date,newdata$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #mario 

#r1
r1start<-as.POSIXct("2011-08-01 00:00:00 CEST")
r1end<-as.POSIXct("2012-03-31 00:00:00 CEST")
#r2
r2start<-as.POSIXct("2012-09-01 00:00:00 CEST")
r2end<-as.POSIXct("2013-02-27 00:00:00 CEST")
#r3
r3start<-as.POSIXct("2013-09-01 00:00:00 CEST")
r3end<-as.POSIXct("2014-02-14 00:00:00 CEST")
#r4
r4start<-as.POSIXct("2014-02-15 00:00:00 CEST")
r4end<-as.POSIXct("2014-08-31 00:00:00 CEST")
#r5
r5start<-as.POSIXct("2014-09-01 00:00:00 CEST")
r5end<-as.POSIXct("2015-02-12 00:00:00 CEST")
#r6
r6start<-as.POSIXct("2015-02-13 00:00:00 CEST")#not 15th feb but 13th feb, and contains about 292 NAs
r6end<-as.POSIXct("2015-09-30 00:00:00 CEST") #changed to sep 30 from aug 31

da[2]<da[1]
classround<-function(date){
  if (is.na(date)) return (NA)
  if (r1start<date & r1end> date) return("r1")
  if (r2start<date & r2end> date) return("r2")
  if (r3start<date & r3end> date) return("r3")
  if (r4start<date & r4end> date) return("r4")
  if (r5start<date & r5end> date) return("r5")
  if (r6start<date & r6end> date) return("r6")
  else return (NA)
}

#da_class<-lapply(da,classround) #not working 
#foo to classify da and fill in a vector 
dclass<-vector()
for (i in 1:length(da)){
  print(i)
  dclass[i]<-classround(da[i])
}
daclass_f<-as.factor(dclass) #curse class as.factor 
summary(daclass_f) #summary 

#summary(daclass_f[newdata$idcol=="E049"]) # a test 
for( i in 1:length(unique(newdata$idcol))){
  id<-unique(newdata$idcol)[i]
  print(id)
  print(summary(daclass_f[newdata$idcol==id]))
}

#######################################################################################################################
######################################################################################################################
#select an ID of GPS collar that at least recorded sufficient information for an animal group in every season.
#create a new variable called Rounds in the data frame to identify relocation points based on the UTC_Date  

newdata$Rounds <- daclass_f
head(newdata)
table(newdata$idcol, newdata$Rounds) #r5, E107 belongs to a different hhid.. confirm whether this claim is justified  

#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E037' & newdata$Rounds=='r1'),] #obs.
newdata.r2 <- newdata[which(newdata$idcol=='E038' & newdata$Rounds=='r2'),]#
newdata.r3 <- newdata[which(newdata$idcol=='e043' & newdata$Rounds=='r3'),]# 
newdata.r4 <- newdata[which(newdata$idcol=='e043' & newdata$Rounds=='r4'),]#
newdata.r5 <- newdata[which(newdata$idcol=='e043' & newdata$Rounds=='r5'),]#
newdata.r6 <- newdata[which(newdata$idcol=='e038' & newdata$Rounds=='r6'),]# 

##special cleaning of round 1

# newdata.r1.1 <- subset(newdata.r1, newdata.r1$Longitude>37.6 & newdata.r1$Longitude<37.75 & newdata.r1$Latitude>3.5 & newdata.r1$Latitude<5.0)
# plot(newdata.r1.1$Longitude, newdata.r1.1$Latitude, col=newdata.r1.1$Rounds)
# newdata.r1 <- newdata.r1.1

hhid.046.herd.f <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.046.herd.f)
plot(hhid.046.herd.f$Longitude, hhid.046.herd.f$Latitude, col=hhid.046.herd.f$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  take lon greater than 37.6 for all households in sarite except for guyo shiriti
hhid.046.herd.w.sub<- subset(hhid.046.herd.w, hhid.046.herd.w$Longitude>38.0 & hhid.046.herd.w$Longitude<39.4 & hhid.046.herd.w$Latitude>4.0 & hhid.046.herd.w$Latitude<4.8)
plot(hhid.046.herd.w.sub$Longitude, hhid.046.herd.w.sub$Latitude, col=hhid.046.herd.w.sub$Rounds)

# cut at lon 39.1 for all the herd in arero

#create a uniquely new ID for these collars.  
hhid.046.herd.w.sub$idcol_new <- as.factor("c013") #add a suffix to as.factor in this line to identify a multiple herds 
hhid.046.herd.f$idcol_new <- as.factor("c013.w") #wara subherd 

#ensure equal length for key variables 
head(hhid.046.herd.w.sub)
length(hhid.046.herd.w.sub$idcol_new) 
length(hhid.046.herd.w.sub$Longitude)
length(hhid.046.herd.w.sub$Latitude) 
length(hhid.046.herd.w.sub$UTC_Date) 
length(hhid.046.herd.w.sub$UTC_Time) 



saveRDS(hhid.046.herd.f, file = "c013.w.hhid.034.wora.subherd.rds")#? fora or wara herd, give a name to specify a herd. 
newdata.34 <- readRDS("c013.w.hhid.034.wora.subherd.rds")
plot(newdata.34$Longitude, newdata.34$Latitude, col=newdata.34$Rounds)

#conclusion: confirm from deployment report if the households 034 and 046 continued herding animals jointly. drop one of them if they herded animals jointly. 
#drop 056 anyway I have a problem with hhid
#################################################################################################################
#################################################################################################################

###hhid 56, name Diqa Kote Abakiyo 

################################################################################################################
################################################################################################################
#import csv files 

list.files()
list_files <- list.files()
length <- length(list_files)

##write cattleTracks to R data frame data type
###identiy a group and animal
start <- proc.time()
newdata<-data.frame()# pre-allocate some variable space for our output data
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

head(newdata)
saveRDS(newdata,file="p.056.allrounds.data.raw.rds") #check if this trajectory is the same as 034 
n.data <- readRDS("p.056.allrounds.data.raw.rds")

sum(is.na(newdata$Longitude)) #a total of 0 nas
sum(is.na(newdata$Latitude)) #a total of 0 nas 
sum(is.na(newdata$Elevation)) #a total of 0 nas
sum(is.na(newdata$UTC_Date))
sum(is.na(newdata$UTC_Time))

# #choose to run these lines only to select complete data without nas

n.data.2 <- n.data[!is.na(n.data$Longitude),]
head(n.data.2)
sum(is.na(n.data.2$Longitude)) #a total of 0 nas
sum(is.na(n.data.2$Latitude)) #a total of 0 nas
sum(is.na(n.data.2$Elevation))
sum(is.na(n.data.2$UTC_Date))
sum(is.na(n.data.2$UTC_Time))
# 
# #new_DF <- DF[is.na(DF$Var),] #an example
newdata <- n.data.2


plot(newdata$Longitude, newdata$Latitude, col=newdata$idcol) #pastoralist had both temporary and permanent place since the inception of study.

da<-strptime(paste(newdata$UTC_Date,newdata$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #mario 

#r1
r1start<-as.POSIXct("2011-08-01 00:00:00 CEST")
r1end<-as.POSIXct("2012-03-31 00:00:00 CEST")
#r2
r2start<-as.POSIXct("2012-09-01 00:00:00 CEST")
r2end<-as.POSIXct("2013-02-27 00:00:00 CEST")
#r3
r3start<-as.POSIXct("2013-09-01 00:00:00 CEST")
r3end<-as.POSIXct("2014-02-14 00:00:00 CEST")
#r4
r4start<-as.POSIXct("2014-02-15 00:00:00 CEST")
r4end<-as.POSIXct("2014-08-31 00:00:00 CEST")
#r5
r5start<-as.POSIXct("2014-09-01 00:00:00 CEST")
r5end<-as.POSIXct("2015-02-12 00:00:00 CEST")
#r6
r6start<-as.POSIXct("2015-02-13 00:00:00 CEST")#
r6end<-as.POSIXct("2015-09-30 00:00:00 CEST") #

da[2]<da[1]
classround<-function(date){
  if (is.na(date)) return (NA)
  if (r1start<date & r1end> date) return("r1")
  if (r2start<date & r2end> date) return("r2")
  if (r3start<date & r3end> date) return("r3")
  if (r4start<date & r4end> date) return("r4")
  if (r5start<date & r5end> date) return("r5")
  if (r6start<date & r6end> date) return("r6")
  else return (NA)
}

#da_class<-lapply(da,classround) #not working 
#foo to classify da and fill in a vector 
dclass<-vector()
for (i in 1:length(da)){
  print(i)
  dclass[i]<-classround(da[i])
}
daclass_f<-as.factor(dclass) #curse class as.factor 
summary(daclass_f) #summary 

#summary(daclass_f[newdata$idcol=="E049"]) # a test 
for( i in 1:length(unique(newdata$idcol))){
  id<-unique(newdata$idcol)[i]
  print(id)
  print(summary(daclass_f[newdata$idcol==id]))
}

#######################################################################################################################
######################################################################################################################
#select an ID of GPS collar that at least recorded sufficient information for an animal group in every season.
#create a new variable called Rounds in the data frame to identify relocation points based on the UTC_Date  

newdata$Rounds <- daclass_f
head(newdata)
table(newdata$idcol, newdata$Rounds) #
#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E042' & newdata$Rounds=='r1'),] #obs.
newdata.r2 <- newdata[which(newdata$idcol=='E042' & newdata$Rounds=='r2'),]#
newdata.r3 <- newdata[which(newdata$idcol=='r064' & newdata$Rounds=='r3'),]# 
newdata.r4 <- newdata[which(newdata$idcol=='e041' & newdata$Rounds=='r4'),]#
newdata.r5 <- newdata[which(newdata$idcol=='e041' & newdata$Rounds=='r5'),]#
newdata.r6 <- newdata[which(newdata$idcol=='e041' & newdata$Rounds=='r6'),]# 



hhid.056.herd.f <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.056.herd.f)
plot(hhid.056.herd.f$Longitude, hhid.056.herd.f$Latitude, col=hhid.056.herd.f$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  take lon greater than 37.6 for all households in sarite except for guyo shiriti
hhid.056.herd.f.sub<- subset(hhid.056.herd.f, hhid.056.herd.f$Longitude>38.96 & hhid.056.herd.f$Longitude<39.35 & hhid.056.herd.f$Latitude>4.0 & hhid.056.herd.f$Latitude<4.7)
plot(hhid.056.herd.f.sub$Longitude, hhid.056.herd.f.sub$Latitude, col=hhid.056.herd.f.sub$Rounds)

# cut at lon 39.1 for all the herd in arero

#create a uniquely new ID for these collars.  
hhid.056.herd.f.sub$idcol_new <- as.factor("c014.f") #add a suffix to as.factor in this line to identify a multiple herds 


#ensure equal length for key variables 
head(hhid.056.herd.f.sub)
length(hhid.056.herd.f.sub$idcol_new) 
length(hhid.056.herd.f.sub$Longitude)
length(hhid.056.herd.f.sub$Latitude) 
length(hhid.056.herd.f.sub$UTC_Date) 
length(hhid.056.herd.f.sub$UTC_Time) 



saveRDS(hhid.056.herd.f.sub, file = "c014.f.hhid.056.fora.subherd.rds")#? fora or wara herd, give a name to specify a herd. 

##create a data for a wara herd 

#create a dataframe for a subherd. hhid 056 

newdata.r1 <-newdata[which(newdata$idcol=='E041' & newdata$Rounds=='r1'),] #obs.
newdata.r2 <- newdata[which(newdata$idcol=='E041' & newdata$Rounds=='r2'),]#
newdata.r3 <- newdata[which(newdata$idcol=='e042' & newdata$Rounds=='r3'),]# 
newdata.r4 <- newdata[which(newdata$idcol=='e042' & newdata$Rounds=='r4'),]#
newdata.r5 <- newdata[which(newdata$idcol=='e042' & newdata$Rounds=='r5'),]#
newdata.r6 <- newdata[which(newdata$idcol=='e042' & newdata$Rounds=='r6'),]# 



hhid.056.herd.w <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.056.herd.w)
plot(hhid.056.herd.w$Longitude, hhid.056.herd.w$Latitude, col=hhid.056.herd.w$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  take lon greater than 37.6 for all households in sarite except for guyo shiriti
hhid.056.herd.w.sub<- subset(hhid.056.herd.w, hhid.056.herd.w$Longitude>38.96 & hhid.056.herd.w$Longitude<39.35 & hhid.056.herd.w$Latitude>4.0 & hhid.056.herd.w$Latitude<4.7)
plot(hhid.056.herd.w.sub$Longitude, hhid.056.herd.w.sub$Latitude, col=hhid.056.herd.w.sub$Rounds)

# cut at lon 39.1 for all the herd in arero

#create a uniquely new ID for these collars.  
hhid.056.herd.w.sub$idcol_new <- as.factor("c014") #add a suffix to as.factor in this line to identify a multiple herds 
summary(hhid.056.herd.w.sub$idcol_new)

#ensure equal length for key variables 
head(hhid.056.herd.w.sub)
length(hhid.056.herd.w.sub$idcol_new) 
length(hhid.056.herd.w.sub$Longitude)
length(hhid.056.herd.w.sub$Latitude) 
length(hhid.056.herd.w.sub$UTC_Date) 
length(hhid.056.herd.w.sub$UTC_Time) 



saveRDS(hhid.056.herd.w.sub, file = "c014.w.hhid.056.wora.subherd.rds")#? fora or wara herd, give a name to specify a herd. 

#################################################################################################################
#################################################################################################################

###hhid 65, name Sasure Shunu  

################################################################################################################
################################################################################################################
#import csv files 

list.files()
list_files <- list.files()
length <- length(list_files)

##write cattleTracks to R data frame data type
###identiy a group and animal
start <- proc.time()
newdata<-data.frame()# pre-allocate some variable space for our output data
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

head(newdata)
#p.730 <- newdata
saveRDS(newdata,file="p.065.allrounds.data.raw.rds") #same trajectory data as 034, a suprise! 

sum(is.na(newdata$Longitude)) #a total of 0 nas
sum(is.na(newdata$Latitude)) #a total of 0 nas 
sum(is.na(newdata$Elevation)) #a total of 0 nas
sum(is.na(newdata$UTC_Date))
sum(is.na(newdata$UTC_Time))

# #choose to run these lines only to select complete data without nas
# n.data <- readRDS("p.540.allrounds.data.raw.rds")
# n.data.2 <- n.data[!is.na(n.data$Longitude),]
# head(n.data.2)
# sum(is.na(n.data.2$Longitude)) #a total of 0 nas
# sum(is.na(n.data.2$Latitude)) #a total of 0 nas
# sum(is.na(n.data.2$Elevation))
# sum(is.na(n.data.2$UTC_Date))
# sum(is.na(n.data.2$UTC_Time))
# 
# #new_DF <- DF[is.na(DF$Var),] #an example
# newdata <- n.data.2


plot(newdata$Longitude, newdata$Latitude, col=newdata$idcol) #pastoralist had both temporary and permanent place since the inception of study.

da<-strptime(paste(newdata$UTC_Date,newdata$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #mario 

#r1
r1start<-as.POSIXct("2011-08-01 00:00:00 CEST")
r1end<-as.POSIXct("2012-03-31 00:00:00 CEST")
#r2
r2start<-as.POSIXct("2012-09-01 00:00:00 CEST")
r2end<-as.POSIXct("2013-02-27 00:00:00 CEST")
#r3
r3start<-as.POSIXct("2013-09-01 00:00:00 CEST")
r3end<-as.POSIXct("2014-02-14 00:00:00 CEST")
#r4
r4start<-as.POSIXct("2014-02-15 00:00:00 CEST")
r4end<-as.POSIXct("2014-08-31 00:00:00 CEST")
#r5
r5start<-as.POSIXct("2014-09-01 00:00:00 CEST")
r5end<-as.POSIXct("2015-02-12 00:00:00 CEST")
#r6
r6start<-as.POSIXct("2015-02-13 00:00:00 CEST")#not 15th feb but 13th feb, and contains about 292 NAs
r6end<-as.POSIXct("2015-09-30 00:00:00 CEST") #changed to sep 30 from aug 31

da[2]<da[1]
classround<-function(date){
  if (is.na(date)) return (NA)
  if (r1start<date & r1end> date) return("r1")
  if (r2start<date & r2end> date) return("r2")
  if (r3start<date & r3end> date) return("r3")
  if (r4start<date & r4end> date) return("r4")
  if (r5start<date & r5end> date) return("r5")
  if (r6start<date & r6end> date) return("r6")
  else return (NA)
}

#da_class<-lapply(da,classround) #not working 
#foo to classify da and fill in a vector 
dclass<-vector()
for (i in 1:length(da)){
  print(i)
  dclass[i]<-classround(da[i])
}
daclass_f<-as.factor(dclass) #curse class as.factor 
summary(daclass_f) #summary 

#summary(daclass_f[newdata$idcol=="E049"]) # a test 
for( i in 1:length(unique(newdata$idcol))){
  id<-unique(newdata$idcol)[i]
  print(id)
  print(summary(daclass_f[newdata$idcol==id]))
}

#######################################################################################################################
######################################################################################################################
#select an ID of GPS collar that at least recorded sufficient information for an animal group in every season.
#create a new variable called Rounds in the data frame to identify relocation points based on the UTC_Date  

newdata$Rounds <- daclass_f
head(newdata)
table(newdata$idcol, newdata$Rounds) #r5, E107 belongs to a different hhid.. confirm whether this claim is justified  

#create a dataframe for a subherd. 

newdata.r1 <-newdata[which(newdata$idcol=='E047' & newdata$Rounds=='r1'),] #obs.
newdata.r2 <- newdata[which(newdata$idcol=='E046' & newdata$Rounds=='r2'),]#
newdata.r3 <- newdata[which(newdata$idcol=='e047' & newdata$Rounds=='r3'),]# 
newdata.r4 <- newdata[which(newdata$idcol=='e047' & newdata$Rounds=='r4'),]#
newdata.r5 <- newdata[which(newdata$idcol=='e047' & newdata$Rounds=='r5'),]#
newdata.r6 <- newdata[which(newdata$idcol=='e046' & newdata$Rounds=='r6'),]# 

##special cleaning of round 1

# newdata.r1.1 <- subset(newdata.r1, newdata.r1$Longitude>37.6 & newdata.r1$Longitude<37.75 & newdata.r1$Latitude>3.5 & newdata.r1$Latitude<5.0)
# plot(newdata.r1.1$Longitude, newdata.r1.1$Latitude, col=newdata.r1.1$Rounds)
# newdata.r1 <- newdata.r1.1

hhid.065.herd.w <- rbind(newdata.r1, newdata.r2, newdata.r3, newdata.r4, newdata.r5, newdata.r6)
head(hhid.065.herd.w)
plot(hhid.065.herd.w$Longitude, hhid.065.herd.w$Latitude, col=hhid.065.herd.w$Rounds)

#clean up newdata for outliers such as obvious mislocation of points due to GPS malfunctioning  take lon greater than 37.6 for all households in sarite except for guyo shiriti
hhid.065.herd.w.sub<- subset(hhid.065.herd.w, hhid.065.herd.w$Longitude>38.9 & hhid.065.herd.w$Longitude<39.3 & hhid.065.herd.w$Latitude>4.0 & hhid.065.herd.w$Latitude<4.8)
plot(hhid.065.herd.w.sub$Longitude, hhid.065.herd.w.sub$Latitude, col=hhid.065.herd.w.sub$Rounds)

# cut at lon 39.1 for all the herd in arero

#create a uniquely new ID for these collars.  
hhid.065.herd.w.sub$idcol_new <- as.factor("c015") #add a suffix to as.factor in this line to identify a multiple herds 
#hhid.056.herd.f$idcol_new <- as.factor("c013.w") #wara subherd 

#ensure equal length for key variables 
head(hhid.065.herd.w.sub)
length(hhid.065.herd.w.sub$idcol_new) 
length(hhid.065.herd.w.sub$Longitude)
length(hhid.065.herd.w.sub$Latitude) 
length(hhid.065.herd.w.sub$UTC_Date) 
length(hhid.065.herd.w.sub$UTC_Time) 



saveRDS(hhid.065.herd.w.sub, file = "c015.hhid.065.fora.subherd.rds")#? fora or wara herd, give a name to specify a herd. 

#####################################################################################################################
####################################################################################################################

### import all data sets -all cattle tracks  

#################################################################################################################
#################################################################################################################

#setwd
list.files()
list_files <- list.files()
CattleTracks <- data.frame()

for(i in 1:length(list_files)){
  print(list_files[i])
  data <- readRDS(list_files[i]) 
  CattleTracks <- rbind(CattleTracks,data)
}

head(CattleTracks[1:20,])
summary(CattleTracks$idcol_new)
#summarysheet <- table(CattleTracks$idcol_new,CattleTracks$Rounds)

saveRDS(CattleTracks, file = "all.cattle.tracks.cleanedV1.rds") #selected and cleaned cattle tracks version1

####################################################################################################################
###################################################################################################################
#### write trajectory object
###################################################################################################################



d2 <- readRDS("all.cattle.tracks.cleanedV1.rds")

sum(is.na(d2$Longitude))
sum(is.na(d2$Latitude))
sum(is.na(d2$Elevation))
length(d2$Longitude)

require(rgdal)
require(raster)


#writing and calling a function
c(d2$Longitude[1],d2$Latitude[1])
getutm<-function(lon,lat){
  xy <- cbind(lon, lat)
  # Convert it to UTM coordinates (in units of meters)
  return(project(xy, "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))#projected coordinate systems for distance calculations.
}

getutm(d2$Longitude[1],d2$Latitude[1]) #pass a matrix into the function

###### begin here to check for duplicated dates in every burst or in a newly created animal id (idcol_new).  
da<-strptime(paste(d2$UTC_Date,d2$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #select dates in the original dataset.
dapos <- as.POSIXct(da, tz="UTC")
da[1]
dapos[1]

animal_date<-paste(as.character(d2$idcol_new),as.character(da)) #paste side by side uniquely id and its dates. 
summary(duplicated(animal_date))## check in the animal_date object if each row is a repeat. yes! returned 216 duplicated entries   
data_wodupes<-d2[-which(duplicated(animal_date)),]#select data without duplicated entries only. 


#confirmation test whether the check for duplicates has been correctly executed. 
length(data_wodupes$Longitude)+216 #should give us the correct length for a column longitude in the original dataset d2. 
length(d2$Longitude)-length(data_wodupes$Longitude)#should give us the number equal to the number of duplicated entries 


locs <-as.data.frame(getutm(d2$Longitude,d2$Latitude)) #select location information in the original dataset. 
head(locs)
tail(locs)
sum(is.na(locs$lon))#na present 0

#drop duplicated entries in da and locs and clean up them to save data in a trajectory object. 

dupes<-which(duplicated(animal_date))#an object containing indexes of duplicated entries (dupes). 
data_wodupes<-d2[-dupes,]            #data without duplicated entries. 

locswithoutdupes<-locs[-dupes,]      #locs without duplicated entries. 
dawithoutdupes<-dapos[-dupes]        #date without duplicated entries. 

#let's create a Spatial Data Frame in UTM zone 87 adding ID, time diff, burst to xy coordinates. 
#creates a spatial data frame from data_wodupes 
#locs <-as.data.frame(getutm(data_wodupes$Longitude,data_wodupes$Latitude)) 

class(locswithoutdupes)
head(locswithoutdupes)
length(locswithoutdupes$lon)

xylocs <- SpatialPoints(locswithoutdupes) #convert to spatial point object, note that it is not a spDF object; then, convert xylocs to spDF
class(xylocs)
proj4string(xylocs) <- CRS("+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

sp.xylocs<- data.frame(xylocs)#create a spatial dataframe from data_wodupes 
idsp <- data.frame(data_wodupes[25])#create a spatial data frame of ID 
dasp <- data.frame(dawithoutdupes)#create a spatial data frame of dawithoutdupes

mergesp <- data.frame(sp.xylocs, idsp, dasp) #merges idcol_new and dawithoutdupes 
coordinates(mergesp) <- 1:2 #convert to spDF object; adds idcol_new and dawithoutdupes  with locations data frame 
plot(mergesp)                #visualize 

speed <- readRDS("all.speed.data_wodupes.rds")
d <- cbind(data_wodupes, speed)

#check out on the lengths -all should be of the same length

length(mergesp$idcol_new) 
length(mergesp$dawithoutdupes) 
length(d$speed)

# length(dasp$dawithoutdupes) 
# length(mergesp$idcol_new) 
# length(mergesp$dawithoutdupes)

require(adehabitatLT)
tr.m1 <- as.ltraj(coordinates(mergesp), date = mergesp$dawithoutdupes, id=mergesp$idcol_new, infolocs = d[,c(1:2, 5:6, 22:24,26)])# 

saveRDS(tr.m1, file = "all.trajectory.data.rds")

head(tr.m1[[1]])
tail(tr.m1[[1]])

###graphical display of the bursts 
plot(tr.m1[1]) #730.bidu wario
plot(tr.m1[2]) #774.Dika Garbicha 
plot(tr.m1[3]) #807.Gufu Galgallo duba fora herd 
plot(tr.m1[4]) #807.Gufu Galgallo duba wara herd 
plot(tr.m1[5]) # 139.nura guyo ilme fora herd
plot(tr.m1[6]) #139.nura guyo ilme wara herd 
plot(tr.m1[7]) #hhid.155.roba kulu fora.subherd 
plot(tr.m1[8]) #hhid.184.dida duba sora .fora.subherd 
plot(tr.m1[9]) #hhid.206.kane karayu. fora.subherd
plot(tr.m1[10]) #hhid.468.jilgaa raqo. wora.subherd 
plot(tr.m1[11]) #hhid.540.jarso qanchora. wora.subherd
plot(tr.m1[12]) #hhid.566.elema galma. fora.subherd
plot(tr.m1[13]) #hhid.566.elema galma. wora.subherd
plot(tr.m1[14]) #hhid.616.guyo shiriti. fora.subherd
plot(tr.m1[15]) #hhid.616.guyo shiriti. wora.subherd 
plot(tr.m1[16]) #hhid.034.tache kote abakiyo. fora.subherd 
plot(tr.m1[17]) #hhid.034.tache kote abakiyo. wora.subherd
plot(tr.m1[18]) #hhid.046.wariyo dido dalacha. fora.subherd
plot(tr.m1[19]) #hhid.046.wariyo dido dalacha. wora.subherd 
plot(tr.m1[20]) #hhid.056.diqa kote. fora.subherd 
plot(tr.m1[21]) #hhid.056.diqa kote. wora.subherd 
plot(tr.m1[22]) #hhid.065.sasure shunu.fora.subherd 



hist(tr.m1[1], "dt", freq=TRUE)
hist(tr.m1[1], "dist", freq=TRUE)

#cutting a burst into several segments on the duration of deployment of GPS collars 

#have a closer look at the values of dt in days according to the dates
plotltr(tr.m1[1], "dt/3600/24")
plotltr(tr.m1[2], "dt/3600/24") 
plotltr(tr.m1[3], "dt/3600/24")
windows()
plotltr(tr.m1[14], "dt/3600/24")
plotltr(tr.m1[15], "dt/3600/24") #a small gps stop has been encountered in april 2014 n restarted after 9 days, data belongs to c011

plotltr(tr.m1[16], "dt/3600/24")#a very good type easily separated into 6 bursts, data belongs to c012
plotltr(tr.m1[17], "dt/3600/24")#can show easily delineable time spans between sucessive GPS deployments 

plotltr(tr.m1[20], "dt/3600/24")
plotltr(tr.m1[21], "dt/3600/24") 
plotltr(tr.m1[22], "dt/3600/24") 

#take a 10-day intervals to segment trajectory into separate bursts. possible to segment trajectory into its respective rounds 

foo <- function(dt) {
  return(dt>(10*3600*24))
}

tr.m1.seg <- cutltraj(tr.m1, "foo(dt)", nextr = TRUE)

#thresholding 

#identify the bursts where the distance between successive relocations was less than 2000 meters at least once. 
bu <- which.ltraj(tr.m1, "dist<2000")
tr.m1[burst(tr.m1)%in%bu$burst] #extract burst satisfying criterion 

##regularize trajectory tr.m1
###regularize a trajectory by placing missing values 
min(mergesp$dawithoutdupes)#returns the oldest date in the data."2011-08-20 06:03:53 UTC"
max(mergesp$dawithoutdupes) #"2015-08-27 16:47:37 UTC" 


refda <- strptime("2011-08-20 00:00:00", "%Y-%m-%d %H:%M:%S") #define a reference date and replace missing values with na  
tr.m1.2 <- setNA(tr.m1, refda, 5, units = "min")

head(tr.m1.2[[2]])
tail(tr.m1.2[[1]])
is.regular(tr.m1.2)#FALSE

refda <- strptime("00:00:00", "%H:%M:%S", tz="UTC") #now set time step to the actual time step
tr.m1.3 <- sett0(tr.m1.2, refda, 5, units = "min") #trajectory object with its time rounded to 5 minutes 

saveRDS(tr.m1.3, file = "all.trajectory.data.RoundedTime.rds")

plotltr(tr.m1.3[1], "dt/3600/24")
plotltr(tr.m1.3[2], "dt/3600/24")

is.regular(tr.m1.3)#TRUE


###data filtering  
################################################################################################################

tr.m1.DF <- ld(tr.m1)  # a data frame of trajectory object

summary(tr.m1.DF$id)

#data cleaning   
tr.m1.DF$dist_km <- tr.m1.DF$dist/1000
tr.m1.DF$dt_hr <- tr.m1.DF$dt/3600
tr.m1.DF$speed_kmh <- tr.m1.DF$dist_km/tr.m1.DF$dt_hr


tr.m1.DF2 <- subset(tr.m1.DF, tr.m1.DF$dt_hr<1 & tr.m1.DF$dist_km<4.5)
length(tr.m1.DF$x)-length(tr.m1.DF2$x)


plot(tr.m1.DF2$dist_km, tr.m1.DF2$dt_hr)
plot(tr.m1.DF2$speed_kmh, tr.m1.DF2$dist_km)

boxplot(tr.m1.DF2$speed_kmh)
boxplot(tr.m1.DF$dt) 
boxplot(tr.m1.DF$dt_hr)
boxplot(tr.m1.DF2$dist_km)


####################################################################################################
##############cut continous variables (speed_kmh) into categorical variables based on velocity categories
####################################################################################################

#testing of chuan's suggest here.
#use cut function, specify bounderies and the resulting values. 
tr.m1.DF2$copyofspeed_kmh <- cut(tr.m1.DF2$speed_kmh, 
                           breaks <- c(-Inf,0.41,1.06,1.94, 2.77, Inf), 
                           labels=c("Stationary", "Heavy Grazing", "Medium Grazing", "Light Grazing", "Travelling"))

table(tr.m1.DF2$copyofspeed_kmh)





#select location information for every sub herd using a sub herd id 
c001 <- subset(tr.m1.DF, id=="c001")
new.c001 <- c001[, c(13, 14)]
write.csv(new.c001, file = "c001.locsdata.csv")


c002 <- subset(tr.m1.DF, id=="c002")
new.c002 <- c002[, c(13, 14)]
head(new.c002)
length(new.c002$Longitude)
write.csv(new.c002, file = "c002.locsdata.csv")

c003 <- subset(tr.m1.DF, id=="c003")
new.c003 <- c003[, c(13, 14)]
head(new.c003)
length(new.c003$Longitude)
write.csv(new.c003, file = "c003.locsdata.csv")

c003.w <- subset(tr.m1.DF, id=="c003.w")
new.c003.w <- c003.w[, c(13, 14)]
head(new.c003.w)
length(new.c003.w$Longitude)
write.csv(new.c003.w, file = "c003_w_locsdata.csv") 

c004 <- subset(tr.m1.DF, id=="c004")
new.c004 <- c004[, c(13, 14)]
head(new.c004)
length(new.c004$Longitude)
write.csv(new.c004, file = "c004.locsdata.csv")


# boxplot(tr.m1.DF$dt)
# plot(tr.m1.DF$dist, tr.m1.DF$dt) # hard to accept travel-distance relationship include a very long distances in kms travelled within a few minutes  
# tr.m1.DF2 <- subset(tr.m1.DF, tr.m1.DF$dist<2000)# dt>260 & dist<1000
# length(tr.m1.DF$x)-length(tr.m1.DF2$x)
# 
# plot(tr.m1.DF2$dist, tr.m1.DF2$dt)
# 
# 
# boxplot(tr.m4DF$dist)
# plot(tr.m4DF$dist, tr.m4DF$rel.angle)
# plot(tr.m4DF$dist, tr.m4DF$speed) 
# 
# 
# 
# ##thresholding of displacement distances  
# thresh<-2000
# subsetbelowthreshold<-tr.m4DF[tr.m4DF$dist<thresh,] #about 652 entries are dropped 
# 
# length(tr.m4DF[,1])-length(subsetbelowthreshold[,1]) #confirm the changes made on the length of our data
# length(subsetbelowthreshold[,1]) +652
# #plot(subsetbelowthreshold$x, subsetbelowthreshold$y)
# 
# tr.m4DF <- subsetbelowthreshold
# tr.m4.2 <- dl(tr.m4DF) #cleaned trajectory. in this trajectory, we maintain all points with a dist <2000m 


#we want to study the trajectory of the day at the scale of the day. we define one trajectory per day. get a daily burst by cutting trajectory
start <- proc.time()
foo <- function(date) {
  da <- as.POSIXlt(date, "UTC")
  ho <- da$hour + da$min/60 
  return(ho>03&ho<03.16)
}

daily.Trips <- cutltraj(tr.m1.3, "foo(date)", nextr = TRUE)

saveRDS(daily.Trips, file = "all.trajectory.data.daily.trips.rds")
proc.time()-start 

is.regular(daily.Trips)

daily.Trips[[1]]
plot(daily.Trips[1]) 
head(daily.Trips[1]) 


daily.Trips[[2]]
plot(daily.Trips[2]) 
head(daily.Trips[2]) 
tail(daily.Trips[[2]]) 


#to view your regular trajectory of points with NA's
summary(daily.Trips)
summary(daily.Trips)[1:10,] 
#now calculating NSD for each point

##import data sets
datansd <- readRDS("all.trajectory.datansd.rds")
daily.Trips <- readRDS("all.trajectory.data.daily.trips.rds")
#d1 <- readRDS("all.trajectory.data.daily.trips.rds")

#d2 <- ld(d1)
#saveRDS(d2, file = "all.trajectory.daily.trips.DF.rds")
d2 <- readRDS("all.trajectory.daily.trips.DF.rds")

#code to get a mean distance walked by herd in every burst (e.g. each trip)
library(devtools)
install_github("ujjwalkarn/xda")
library(xda)
#numSummary(iris)
#https://www.r-bloggers.com/introducing-xda-r-package-for-exploratory-data-analysis/ 
#summary.d2 <- numSummary(data.frame(d2[,6])) #

summary.d3 <- aggregate(d2$dist, by=list(d2$burst), na.rm=TRUE, FUN=sum)#return a sum for (dist) by group (burst)
colnames(summary.d3) <- c("burst","Path_dist") 
summary.d4 <- data.frame(summary.d3[,2])
colnames(summary.d4) <- "Path_dist"          # a trip length in every trip gone

trips <- data.frame(summary(daily.Trips))    #convert a summary of ltraj daily trip to dataframe
trips <- data.frame(trips,summary.d4)        #join trip length and the trips 
boxplot(trips$Path_dist) #some outliers. select daily trajectory a cut off point at less than 40km 
trips2[200:280,] # path_dist has got 0 values, a strange thing! yeap! noticed na values 

##delineation of movement periods in trips 
##assign a year as when we would expect our herds to disperse or migrate from wet to dry season grazing areas. 
trips2 <- trips
range(trips$date.begin) #"2011-08-20 06:05:00 UTC" "2015-08-27 03:10:00 UTC" 

trips2$Year <- NULL
trips2$Year[trips2$date.begin>"2011-08-20 06:00:00 UTC" & trips2$date.begin<"2011-09-16 03:10:00 UTC" ] <- "LD2011"
trips2$Year[trips2$date.begin>"2011-09-15 03:10:00 UTC" & trips2$date.begin<"2012-01-01 03:10:00 UTC" ] <- "SR2011"
trips2$Year[trips2$date.begin>"2011-12-31 03:10:00 UTC" & trips2$date.begin<"2012-03-01 03:10:00 UTC" ] <- "SD2012"
trips2$Year[trips2$date.begin>"2012-02-29 03:10:00 UTC" & trips2$date.begin<"2012-06-01 03:10:00 UTC" ] <- "LR2012"
trips2$Year[trips2$date.begin>"2012-05-31 03:10:00 UTC" & trips2$date.begin<"2012-09-16 03:10:00 UTC" ] <- "LD2012"
trips2$Year[trips2$date.begin>"2012-09-15 03:10:00 UTC" & trips2$date.begin<"2013-01-01 03:10:00 UTC" ] <- "SR2012"
trips2$Year[trips2$date.begin>"2012-12-31 03:10:00 UTC" & trips2$date.begin<"2013-03-01 03:10:00 UTC" ] <- "SD2013"
trips2$Year[trips2$date.begin>"2013-02-28 03:10:00 UTC" & trips2$date.begin<"2013-06-01 03:10:00 UTC" ] <- "LR2013"
trips2$Year[trips2$date.begin>"2013-05-31 03:10:00 UTC" & trips2$date.begin<"2013-09-16 03:10:00 UTC" ] <- "LD2013"
trips2$Year[trips2$date.begin>"2013-09-15 03:10:00 UTC" & trips2$date.begin<"2014-01-01 03:10:00 UTC" ] <- "SR2013"
trips2$Year[trips2$date.begin>"2013-12-31 03:10:00 UTC" & trips2$date.begin<"2014-03-01 03:10:00 UTC" ] <- "SD2014"
trips2$Year[trips2$date.begin>"2014-02-28 03:10:00 UTC" & trips2$date.begin<"2014-06-01 03:10:00 UTC" ] <- "LR2014"
trips2$Year[trips2$date.begin>"2014-05-31 03:10:00 UTC" & trips2$date.begin<"2014-09-16 03:10:00 UTC" ] <- "LD2014"
trips2$Year[trips2$date.begin>"2014-09-15 03:10:00 UTC" & trips2$date.begin<"2015-01-01 03:10:00 UTC" ] <- "SR2014"
trips2$Year[trips2$date.begin>"2014-12-31 03:10:00 UTC" & trips2$date.begin<"2015-03-01 03:10:00 UTC" ] <- "SD2015"
trips2$Year[trips2$date.begin>"2015-02-28 03:10:00 UTC" & trips2$date.begin<"2015-06-01 03:10:00 UTC" ] <- "LR2015"
trips2$Year[trips2$date.begin>"2015-05-31 03:10:00 UTC" & trips2$date.begin<"2015-09-16 03:10:00 UTC" ] <- "LD2015"
#trips2$Year[trips2$date.begin>"2015-09-15 03:10:00 UTC" & trips2$date.begin<"2016-01-01 03:10:00 UTC" ] <- SR2015

trips2$Year <- as.factor(trips2$Year)
boxplot(trips2$Path_dist~trips2$Year) #drop all observations with a higher na values for variable NAS. 
saveRDS(trips2, file = "all.trajectory.daily.trips.DF.Tripsummary.rds")

#data preprocessing chain to apply a k-means clustering    

range(d2$date) #"2011-08-20 06:05:00 UTC" "2015-08-27 16:50:00 UTC"

d2$Year <- NULL
d2$Year[d2$date>"2011-08-20 06:00:00 UTC" & d2$date<"2011-09-16 03:10:00 UTC" ] <- "LD2011"
d2$Year[d2$date>"2011-09-15 03:10:00 UTC" & d2$date<"2012-01-01 03:10:00 UTC" ] <- "SR2011"
d2$Year[d2$date>"2011-12-31 03:10:00 UTC" & d2$date<"2012-03-01 03:10:00 UTC" ] <- "SD2012"
d2$Year[d2$date>"2012-02-29 03:10:00 UTC" & d2$date<"2012-06-01 03:10:00 UTC" ] <- "LR2012"
d2$Year[d2$date>"2012-05-31 03:10:00 UTC" & d2$date<"2012-09-16 03:10:00 UTC" ] <- "LD2012"
d2$Year[d2$date>"2012-09-15 03:10:00 UTC" & d2$date<"2013-01-01 03:10:00 UTC" ] <- "SR2012"
d2$Year[d2$date>"2012-12-31 03:10:00 UTC" & d2$date<"2013-03-01 03:10:00 UTC" ] <- "SD2013"
d2$Year[d2$date>"2013-02-28 03:10:00 UTC" & d2$date<"2013-06-01 03:10:00 UTC" ] <- "LR2013"
d2$Year[d2$date>"2013-05-31 03:10:00 UTC" & d2$date<"2013-09-16 03:10:00 UTC" ] <- "LD2013"
d2$Year[d2$date>"2013-09-15 03:10:00 UTC" & d2$date<"2014-01-01 03:10:00 UTC" ] <- "SR2013"
d2$Year[d2$date>"2013-12-31 03:10:00 UTC" & d2$date<"2014-03-01 03:10:00 UTC" ] <- "SD2014"
d2$Year[d2$date>"2014-02-28 03:10:00 UTC" & d2$date<"2014-06-01 03:10:00 UTC" ] <- "LR2014"
d2$Year[d2$date>"2014-05-31 03:10:00 UTC" & d2$date<"2014-09-16 03:10:00 UTC" ] <- "LD2014"
d2$Year[d2$date>"2014-09-15 03:10:00 UTC" & d2$date<"2015-01-01 03:10:00 UTC" ] <- "SR2014"
d2$Year[d2$date>"2014-12-31 03:10:00 UTC" & d2$date<"2015-03-01 03:10:00 UTC" ] <- "SD2015"
d2$Year[d2$date>"2015-02-28 03:10:00 UTC" & d2$date<"2015-06-01 03:10:00 UTC" ] <- "LR2015"
d2$Year[d2$date>"2015-05-31 03:10:00 UTC" & d2$date<"2015-09-16 03:10:00 UTC" ] <- "LD2015"
#trips2$Year[d2$date>"2015-09-15 03:10:00 UTC" & d2$date<"2016-01-01 03:10:00 UTC" ] <- SR2015

d2$Year <- as.factor(d2$Year)
boxplot(d2$dist~d2$Year) # exceptionally longest steps! drop them. 

datansd2 <- data.frame(datansd[,19:24])
d3 <- data.frame(d2, datansd2)

saveRDS(d3, file = "all.trajectory.datansd.final.rds")
d3.2 <- na.omit(d3)                                     #cleaned dataset for analysis. 
boxplot(d3.2$dist~d3.2$Year)

##thresholding 
thresh<-3500
subsetbelowthreshold<-d3.2[d3.2$dist<thresh,]

length(d3.2[,1])-length(subsetbelowthreshold[,1])
length(subsetbelowthreshold[,1]) 
#plot(subsetbelowthreshold$Longitude, subsetbelowthreshold$Latitude)

d3.3 <- subsetbelowthreshold
boxplot(d3.3$dist~d3.3$Year)
saveRDS(d3.3, file = "all.traj.cleaned.thesholdbelow3500.rds")

d4 <- readRDS("all.traj.cleaned.thesholdbelow3500.rds")
plot(d4$minsday, d4$R2n)#daylight activity cut-off 240-800minsday
range(d4$minsday)#0 1435 
range(d4$hdaylitr2)#0 23

d4$minsday2 <- as.factor(d4$hdaylitr2) 
boxplot(d4$dist~d4$minsday2, xlab="Hours in UTC", ylab="Step lengths (m)", main="Step lengths vs. Time of day") 

d4.daylights <- subset(d4, d4$hdaylitr2>3 & d4$hdaylitr2<16)
d4.daylights.outward <- subset(d4, d4$hdaylitr2>3 & d4$hdaylitr2<8)
d4.daylights.res.extraction <- subset(d4, d4$hdaylitr2>7 & d4$hdaylitr2<13)
d4.daylights.return <- subset(d4, d4$hdaylitr2>12 & d4$hdaylitr2<16)

length(d4.daylights.outward$x)+length(d4.daylights.res.extraction$x)+length(d4.daylights.return$x)
range(d4.daylights.outward$hdaylitr2) 
range(d4.daylights.res.extraction$hdaylitr2) 
range(d4.daylights.return$hdaylitr2) 
range(d4.daylights$hdaylitr2) # 4 15

length(d4$x)-length(d4.daylights$x)
length(d4.daylights$x)+ 2116828


plot(d4.daylights.res.extraction$zero1, d4.daylights.res.extraction$R2n) 
boxplot(d4.daylights.res.extraction$dist~d4.daylights.res.extraction$minsday2, xlab="Hours in UTC", ylab="Step lengths (m)", main="Step lengths vs. Time of day") 

##thresholding R2N 
thresh <- 500000000 #about 14km
d4.daylights.res.extraction.1 <- d4.daylights.res.extraction[d4.daylights.res.extraction$R2n<thresh,]
plot(d4.daylights.res.extraction.1$zero1, d4.daylights.res.extraction.1$R2n) 

length(d4.daylights.res.extraction$x)-length(d4.daylights.res.extraction.1$x)#22 observations dropped. 
length(d4.daylights.res.extraction.1$x)+ 22
boxplot(d4.daylights.res.extraction.1$dist~d4.daylights.res.extraction.1$Year)

d4.daylights.res.extraction.2<- subset(d4.daylights.res.extraction.1, select=c("dist","rel.angle","R2n")) #subset 
d4.daylights.res.extraction.2$distlog10 <- log10(d4.daylights.res.extraction.2$dist)
d4.daylights.res.extraction.2$R2nlog10 <- log10(d4.daylights.res.extraction.2$R2n+1)
d4.daylights.res.extraction.3 <- subset(d4.daylights.res.extraction.2, select = c("rel.angle", "distlog10", "R2nlog10")) 

library(PerformanceAnalytics)
chart.Correlation(d4.daylights.res.extraction.3, 
                  method = "pearson",
                  histogram=TRUE,
                  adjust="none",
                  alpha=0.05)

d4.daylights.res.extraction.4 <- scale(d4.daylights.res.extraction.3)         #subset scaled data 

#d4 <- scale(d4[-1]) #drop first column
#d5 <- subset(d3.3, select=c("dist","rel.angle","R2n"))  #subset
                                   
wss <- (nrow(d4.daylights.res.extraction.4))*sum(apply(d4.daylights.res.extraction.4,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(d4.daylights.res.extraction.4,iter.max = 100, algorithm="MacQueen",
                                     centers = i)$withinss)

plot(1:15, wss, type="b", xlab="Number of clusters", ylab="Within group sum of squares")

#kmean cluster analysis
require(cluster)
km <-kmeans(d4.daylights.res.extraction.4,8, iter.max = 50) 
if (km$ifault==4) { km = kmeans(d4.daylights.res.extraction.4,iter.max = 50, km$centers, algorithm="MacQueen"); }

attributes(km)
km$centers
km$size #we can look at the composition of clusters  
km$withinss #We can look at the within sum of squares of each cluster

aggregate(d4.daylights.res.extraction.2, by=list(km2$cluster), FUN=mean)##We can compare summary statistics.

#Group 1. STATIONARY: a very short segments, and at home. may be the animal or collar stayed at home (68152).
#Group 2. HEAVY GRAZING: slow-paced walks, highly tortous, and about 4166 meters away from the starting point (187302). 
#Group 3. HEAVY GRAZING: slow-paced walks, highly tortous, and about 4177 meters away from the starting point (186001).
#Group 4. WALKING: long-paced walks, persistently directional and 5054 away from the starting point (size 443989). include trips to water points and journey to farfields. here, the persistence directions carries more weight.    


#append cluster assignment 
d4.daylights.res.extraction.5 <- data.frame(d4.daylights.res.extraction.1, km$cluster) #the thresholded dataset
d4.daylights.res.extraction.5$km.cluster <- as.factor(d4.daylights.res.extraction.5$km.cluster)
plot(d4.daylights.res.extraction.5$x, d4.daylights.res.extraction.5$y, pch=20, col=d4.daylights.res.extraction.5$km.cluster)
##############################
##split group 4 further into 5 sub-groups
#################################
d4.daylights.res.extraction.4.grp4 <- data.frame(d4.daylights.res.extraction.4, km$cluster) 
d4.daylights.res.extraction.4.grp4 <- d4.daylights.res.extraction.4.grp4[d4.daylights.res.extraction.4.grp4$km.cluster==4,] 
d4.daylights.res.extraction.4.grp4.sub <- subset(d4.daylights.res.extraction.4.grp4, select = c("rel.angle", "distlog10", "R2nlog10"))

wss <- (nrow(d4.daylights.res.extraction.4.grp4.sub))*sum(apply(d4.daylights.res.extraction.4.grp4.sub,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(d4.daylights.res.extraction.4.grp4.sub,iter.max = 100, algorithm="MacQueen",
                                     centers = i)$withinss)
plot(1:15, wss, type="b", xlab="Number of clusters", ylab="Within group sum of squares")

km <-kmeans(d4.daylights.res.extraction.4.grp4.sub,5, iter.max = 50) 
if (km$ifault==4) { km = kmeans(d4.daylights.res.extraction.4.grp4.sub,iter.max = 50, km$centers, algorithm="MacQueen"); }

km$size 

d4.daylights.res.extraction.5.grp4 <- d4.daylights.res.extraction.5[d4.daylights.res.extraction.5$km.cluster==4,]
d4.daylights.res.extraction.5.grp4.sub <- subset(d4.daylights.res.extraction.5.grp4, select = c("dist", "R2n", "rel.angle"))
aggregate(d4.daylights.res.extraction.5.grp4.sub, by=list(km$cluster), FUN=mean)##We can compare summary statistics.




#biological interpretation for daylights locations (ALL)
#Group 1. GRAZING: slow-paced walks, highly tortous, happened at around 2758 meters away from the starting points and occured on average timing of 6.9 hours into the daily grazing itinerary. State 1 is found in matatika. *****

#group 2. slow-paced walks, persistently directional, happened at a distant place about 3003 meters from the starting point, and happened in the early hours of the grazing itinerary at about 3.65 hours. State 2 is an outward trip to the field. may represent a light grazing as the herds approaches the main graing fields. ++++ 

#Group 3. Walking: a relatively fast-paced walks, persistently directional, found in the distant places like 8821 meters away from the starting point, and happened around 6.84 hours since the start of the journey. the timing of this state is within the peak activity period. in matatika 

#Group 4. GRAZING: a slow-paced walks, highly tortous, about 2772 meters away from the starting point, and 6.9 hours hours into grazing itinerary. movement characteristcs of the state 4 are quite similar to those of state 1. *****  

#Group 5. a short-paced walks, persistently directional, at about 2861 meters to the starting point and late  into a daily itinerary at about 9.9 hours. state 5 is part of a return trip from the field and to a starting place. ++++ 

#Group 6. Walking:  a long-paced steps, persistently directional, about 3164.8 meters to a starting place, and way too late in the evening at about 10.6 hours into daily itinerary. state 6 is part of a return trip to a starting place.

#Group 7. walking: a long-paced steps, persistenly directional, about 3328 meters away from the starting point and only 2.9 hours into the morning session. state 7 is part of an outward trip.  

###################################
#thresholding nsd <2e+08
#group 1 fast-paced persistent walks, happened within the range of 3467 meters from the starting point, and about 6.33 hours elapsed since the journey started. the state happened at the peak grazing hours. IN matatika. +++

#group 2 GRAZING, a short-paced walks, highly tortous- occurs at about 2684 meters away from the starting point, at the peak hours about 6.8 hours. In matatika *****

#group 3 slightly slow-paced movements, persitent walking on a return journey approaching the starting point and in the range of of 2700 meters to a starting point, & in the evening at 10 hours. maruu galchuum. XXX

#Group 4 fast-paced persistent walks, happened in the farfield about 8143 meters from the starting point, at the peak hours around 7 hours since the start of the journey. BEYOND MATATIKA: A journey to and from a water point or a movement state completed during migration. +++  

#Group 5 GRAZING a short-paced walks, highly tortous, happened at about 2670 meters from a start place, and about 7 hours since the journey started. In a matatika. movement ccs of the state 5 is quite similar with the state 2 ***** 

#Group 6 slightly slow-paced walks, an outward movement in the closest range of the starting point about 2833 meters, and happened just within the first 3 hours into the start of the journey. Boobaya. XXXX


####################################################################################################################
######################################################################################################################

dem<-raster("D:\\GPS\\A_DEM_mask3.tif")
veg <- raster("D:/shibia/herds/a_vegdevVer9.tif")
extent(veg)
plot(veg)

image(veg)



#
#task:     Extracting cell values from an environmental raster (extract dem values)  
#using a series of corresponding points (e.g. at herd location points) 
#Use:     make use of geocomputaion capabilities of SAGA GIS from within R

# ------------------------------------------------------------
# Initial settings:
# ------------------------------------------------------------
install.packages("RSAGA")
library(RSAGA); rsaga.env()

#myenv <- rsaga.env(workspace="C:/temp/", path="c:/Program Files/QGIS Wien/apps/saga")#
myenv <- rsaga.env(workspace="C:/temp/", path="c:/temp/saga-gis")# 

# rsaga.esri.to.sgrd(in.grids="tmp.asc", 
#                    out.sgrds="Zn_RK.sgrd", in.path=getwd(),env=myenv) 


##cattle herd location data
# setwd("D:/shibia/herds") 
# d1<-readRDS("cattleGPSfiles")#to read a .RDS file

# d1$Longitude<-as.numeric(d1$Longitude)
# d1$Elevation<-as.numeric(d1$Elevation)
# sum(is.na(d1$Longitude))
# 
# d2 <- d1[, c(1:12, 22,23)]#subset
# head(d2)


#d2 <- na.omit(d2, cols=c("Longitude", "Latitude", "Elevation"))## omit rows where either 'x', 'y' or 'z' have missing values
# d2$UTC_Date <- as.numeric(d2$UTC_Date)
# d2$UTC_Date <- as.numeric(d2$UTC_Time)
#d2$groupcol <- as.factor(d2$groupcol)
#d2$idcol <- as.factor(d2$idcol)
is.factor(d1$groupcol)
is.factor(d1$idcol)

require(rgdal)
require(raster)
#d1$Longitude[1],d1$Latitude[1]

#writing and calling a function
c(d2$Longitude[1],d2$Latitude[1])
getutm<-function(lon,lat){
  xy <- cbind(lon, lat)
  # Convert it to UTM coordinates (in units of meters)
  return(project(xy, "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))#projected coordinate systems for distance calculations.
}
getdist<-function(point1,point2){
  dv<-(c(point1[1]-point2[1],point1[2]-point2[2],point1[3]-point2[3]))
  return(sqrt(dv[1]^2+dv[2]^2+dv[3]^2))
}

#pass a matrix into the function
getutm(d2$Longitude[1],d2$Latitude[1])

###################################################################################################################
##################################################################################################################
### create an object of class trajectory
##################################################################################################################
#################################################################################################################
#d1
head(d1)
d1$place <- getutm(d1$Longitude,d1$Latitude) 
#points(d1$place)
plot(d1$place, pch=20, col=d1$groupcol)

###############################################################################################################
##############################################################################################################
#######First the date needs to be transformed into an object of the class POSIXct. 
#############################################################################################################
#############################################################################################################
class(d1$UTC_Date)#factor
#da <- mdy(d1$UTC_Date, tz="UTC")#date conversion 
da<-strptime(paste(d1$UTC_Date,d1$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #mario 
#da <- as.POSIXct(strptime(as.character(d1$UTC_Date),"%m%d%y")) #We  pkg version
da <- as.POSIXct(da, tz="UTC") 
#da <- as.POSIXct(da, tz = "GMT")
#da <- as.Date(d1$UTC_Date)

#da <- as.character(d1$UTC_Date)
head(da)
tail(da)
class(da)
#locs$Date <- as.POSIXct(locs$Date)
#da <- as.POSIXct(strptime(as.character(d1$UTC_Date),"%y%m%d")) #
#d4$UTC_Date <- mdy(d4$UTC_Date, tz="Africa/Nairobi")#date conversion
# #d2$UTC_Date <- as.POSIXct(strptime(as.character(d2$UTC_Date),"%m%d%y"))
#c(d4$Longitude[1],d4$Latitude[1]) 
locs <- d1$place
head(locs)
#We can then create an object of class ltraj to store the herd movements: 
tr1 <- as.ltraj(coordinates(locs), date = da, id = d1$idcol)#function computes the descriptive parameters from X and Y and the date. degrees are expressed in radians 
class(tr1)
head(tr1[[1]])
tail(tr1[[1]])

head(tr1[[2]]) 
##graphical display of the bursts
plot(tr1)#returns an error
tr.round1 <- ld(tr1) #converting trajectory to data frame
head(tr.round1) #an excellent dataframe 

saveRDS(tr.round1,file="tr.round1.rds")  
dl(tr.round1)#converting a data frame to a trajectory type object. 

is.regular(tr1)#false

#write summary to excel file 
summary(tr1)
install.packages("xlsx")
library(xlsx)
write.xlsx(a, "trajectory_round1_summary.xlsx")#write a trajectory summary in a spread sheet.  

####################################################################################################################
###################################################################################################################
###############Cutting a burst into several segments 
##################################################################################################################
##################################################################################################################




################################################################################################################
################################################################################################################
#### subset round 1 data; create a subset herd for every pastoralist in a given season.  
################################################################################################################
################################################################################################################
#Dika Garbicha Jattani: siqu of yabello distric
ID774.1 <-d2[which(d2$groupcol=='774' & d2$idcol=='E052'),] #17717 obs.
ID774.2 <-d2[which(d2$groupcol=='774' & d2$idcol=='E053'),] #27516 obs.
ID774.3 <-d2[which(d2$groupcol=='774' & d2$idcol=='E054'),] #46903 obs. 

#ID774 <- d2[d2$groupcol=="774",] 

head(ID774.3)
d3 <- ID774.3  #at the end, this object is renamed to become d4. 

require(rgdal)
require(raster)
#d1$Longitude[1],d1$Latitude[1]

#writing and calling a function
c(d1$Longitude[1],d1$Latitude[1])
getutm<-function(lon,lat){
  xy <- cbind(lon, lat)
  # Convert it to UTM coordinates (in units of meters)
  return(project(xy, "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))#projected coordinate systems for distance calculations.
}

getutm(d1$Longitude[1],d1$Latitude[1])


dem<-raster("D:\\GPS\\A_DEM_mask3.tif")
image(dem)

d3$place <- getutm(d3$Longitude,d3$Latitude) 
points(d3$place)

#plot(d3$place, type="l")#misplaced gps locations are present in the data.
#points(locs, pch = 20, col = locs$Name) 
#plot(d3$Longitude,d3$Latitude,type="l")



#########################################################################################
####extracting elevation data at herd location points and use these points to calculate speeds of herds.  
####http://neondataskills.org/R/extract-raster-data-R/ (a reference)
####https://anonymousthis.wordpress.com/2013/02/12/extracting-point-values-from-rasters-using-r-gdal-by-example/ 
####http://www.organicdatascience.org/gleonfellowship/index.php/Zonal_Statistics_or_Extracting_Raster_Data 
##########################################################################################
require(rgdal)
require(raster)
#d1$Longitude[1],d1$Latitude[1]

#writing and calling a function
c(d1$Longitude[1],d1$Latitude[1])
getutm<-function(lon,lat){
  xy <- cbind(lon, lat)
  # Convert it to UTM coordinates (in units of meters)
  return(project(xy, "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))#projected coordinate systems for distance calculations.
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
p1<-c(getutm(d3$Longitude[1],d3$Latitude[1]),d3$Elevation[1])
p2<-c(getutm(d3$Longitude[2],d3$Latitude[2]),d3$Elevation[2])
getdist(p1,p2)

# d1$dist<-rep(0,length(d1[,1]))
# for(i in 1:length(d1[,1])){
#   p1<-c(getutm(d1$Longitude[i],d1$Latitude[i]),d1$Elevation[i])
#   p2<-c(getutm(d1$Longitude[i+1],d1$Latitude[i+1]),d1$Elevation[i+1])
#   d1$dist[i+1]<-getdist(p1,p2)
# }
# 
# hist(d1$Elevation)
# plot(d1$Elevation[1:1000],type="l")
# install.packages("signal")
# require(signal)
# 
# plot(d1$Elevation[1:1000],type="l")
# plot(sgolayfilt(d1$Elevation[1:1000],3,7,0,1),type="l")
#dem<-raster("D:\\shibia\\herds\\dem_msk_resampled") 
#dem<-raster("D:\\PAPER1\\A_DEM_mask.tif")
extent(dem)
crs(dem)
plot(dem)
extract(dem,cbind(500000,500000))

start <- proc.time() 
d3$Ele_DEM<-rep(0,length(d3[,1]))
for(i in 1:length(d3[,1])){
  p<-getutm(d3$Longitude[i],d3$Latitude[i])
  d3$Ele_DEM[i]<-extract(dem,cbind(p[1],p[2]))
}

proc.time()-start 
# for(i in 1:100){
#   p<-getutm(d1$Longitude[i],d1$Latitude[i])
#   d1$Ele_DEM[i]<-extract(dem,cbind(p[1],p[2]))
# }

## extract(dem,cbind(500000,500000))
## getutm(d1$Longitude[1],d1$Latitude[1])
## getutm(d1$Longitude[10],d1$Latitude[10])
## d1[10,]


####
p1<-c(getutm(d3$Longitude[1],d3$Latitude[1]),d3$Ele_DEM[1])
p2<-c(getutm(d3$Longitude[2],d3$Latitude[2]),d3$Ele_DEM[2])
getdist(p1,p2)

start <- proc.time() 
d3$dist2<-rep(0,length(d3[,1]))
for(i in 1:100){
  p1<-c(getutm(d3$Longitude[i],d3$Latitude[i]),d3$Elevation[i])
  p2<-c(getutm(d3$Longitude[i+1],d3$Latitude[i+1]),d3$Elevation[i+1])
  d3$dist2[i+1]<-getdist(p1,p2)
}
proc.time()-start

start <- proc.time() 
d3$dist.test<-rep(0,length(d3[,1]))
for(i in 1:46903){
  p1<-c(getutm(d3$Longitude[i],d3$Latitude[i]),d3$Ele_DEM[i])
  p2<-c(getutm(d3$Longitude[i+1],d3$Latitude[i+1]),d3$Ele_DEM[i+1])
  d3$dist.test[i+1]<-getdist(p1,p2)
}
proc.time()-start
######################################################################################################
####time
#######################################################################################################
d3$UTC_Time
paste(d3$UTC_Date[1],d3$UTC_Time[1],sep=" ")
date1<-strptime(paste(d3$UTC_Date[1],d3$UTC_Time[1],sep=" "),format="%m/%d/%Y %H:%M:%S")
date2<-strptime(paste(d3$UTC_Date[2],d3$UTC_Time[2],sep=" "),format="%m/%d/%Y %H:%M:%S")

dt<-difftime(date2,date1,units="secs")
as.numeric(dt)

start <- proc.time()
d3$deltat_seconds<-rep(0,length(d3[,1]))
d3$speed<-rep(0,length(d3[,1]))
for (i in 1:46903){
  date1<-strptime(paste(d3$UTC_Date[i],d3$UTC_Time[i],sep=" "),format="%m/%d/%Y %H:%M:%S")
  date2<-strptime(paste(d3$UTC_Date[i+1],d3$UTC_Time[i+1],sep=" "),format="%m/%d/%Y %H:%M:%S")
  dt<-difftime(date2,date1,units="secs")
  d3$deltat_seconds[i+1]<-as.numeric(dt)
  d3$speed<-d3$dist.test/as.numeric(dt)
}
proc.time()-start
# ?difftime
# 
# as.difftime(c(as.character(d1$UTC_Time[1]),as.character(d1$UTC_Time[2])))

##############################################################################################################
####thresholding 
#############################################################################################################

# par(mfrow=c(2,1))
# plot(d1$Elevation[1:7200],type="l")
# plot(d1$Ele_DEM[1:5000],type="l")
# 600*5
# d1$date<-rep(as.POSIXct(time),length(d1[,1]))
# for(i in 1:length(d1[,1])){
#   time<-strptime(paste(d1$UTC_Date[i],d1$UTC_Time[i],sep=" "),format="%m/%d/%Y %H:%M:%S")
#   d1$date[i]<-as.POSIXct(time)
#   ###d1$date[i]<-as.numeric(strftime(time,format="%H"))
#   
# }
# par(mfrow=c(2,1))
# start<-100
# range<-40000
# ###plot(d1$date[1:range],d1$Elevation[1:range],type="l")
# plot(d1$Ele_DEM[start:(start+range)],type="l")
# plot(d1$date[start:(start+range)],d1$Ele_DEM[start:(start+range)],type="l")
# plot(d1$Longitude,d1$Latitude,type="l")
# hist(d1$dist2[1:500])
# plot(d1$dist2,log="y")

# breaks<-c(seq(0,1000,by=200),seq(2000,100000,by=10000),max(na.exclude(d1$dist2))+1)
# hist(na.exclude(d1$dist2),breaks=500)
# max(na.exclude(d1$dist2))
thresh<-2000
subsetbelowthreshold<-d3[d3$dist.test<thresh,]
length(d3[,1])-length(subsetbelowthreshold[,1])
length(subsetbelowthreshold[,1]) 
plot(subsetbelowthreshold$Longitude, subsetbelowthreshold$Latitude)

saveRDS(subsetbelowthreshold,file="ID774.herd3.rds") 
l <- readRDS("ID774.herd3.rds")
#write.csv(subsetbelowthreshold,file="ID774.herd3")


R1.774.e054 <- readRDS("ID774.herd3.rds")
head(R1.774.e054)
plot(R1.774.e054$Longitude, R1.774.e054$Latitude)
d4 <- R1.774.e054


################################################################################################################
###############################################################################################################
################        data round 2
###############################################################################################################
###############################################################################################################
#import csv files 
list.files()
list_files <- list.files()
length <- length(list_files)

##write cattleTracks to R data frame data type
###identiy a group and animal
start <- proc.time()
newdata<-data.frame()# pre-allocate some variable space for our output data
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

head(newdata)
R2 <- newdata
saveRDS(newdata,file="Round2.rds") 

################################################################################################################
#writing and calling a function to calculate dist, and check for outliers in the relocation data 
###############################################################################################################
R2$place <- getutm(R2$Longitude,R2$Latitude) 
c(R2$Longitude[1],R2$Latitude[1])

#writing and calling a function

getutm<-function(lon,lat){
  xy <- cbind(lon, lat)
  return(project(xy, "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))#projected coordinate systems for distance calculations.
}
getdist<-function(point1,point2){
  dv<-(c(point1[1]-point2[1],point1[2]-point2[2],point1[3]-point2[3]))
  return(sqrt(dv[1]^2+dv[2]^2+dv[3]^2))
}

#pass a matrix into the function
getutm(R2$Longitude[1],R2$Latitude[1])


p1<-c(getutm(R2$Longitude[1],R2$Latitude[1]),R2$Elevation[1])
p2<-c(getutm(R2$Longitude[2],R2$Latitude[2]),R2$Elevation[2])
getdist(p1,p2)

start <- proc.time() 
R2$dist2<-rep(0,length(R2[,1]))
for(i in 1:1842883){
  p1<-c(getutm(R2$Longitude[i],R2$Latitude[i]),R2$Elevation[i])
  p2<-c(getutm(R2$Longitude[i+1],R2$Latitude[i+1]),R2$Elevation[i+1])
  R2$dist2[i+1]<-getdist(p1,p2)
}
proc.time()-start

######################################################################################################
####time
#######################################################################################################

#R2$UTC_Time
paste(R2$UTC_Date[1],R2$UTC_Time[1],sep=" ")
date1<-strptime(paste(R2$UTC_Date[1],R2$UTC_Time[1],sep=" "),format="%m/%d/%Y %H:%M:%S")
date2<-strptime(paste(R2$UTC_Date[2],R2$UTC_Time[2],sep=" "),format="%m/%d/%Y %H:%M:%S")

dt<-difftime(date2,date1,units="secs")
as.numeric(dt)

start <- proc.time()
R2$deltat_seconds<-rep(0,length(R2[,1]))
R2$speed<-rep(0,length(R2[,1]))
for (i in 1:1842883){
  date1<-strptime(paste(R2$UTC_Date[i],R2$UTC_Time[i],sep=" "),format="%m/%d/%Y %H:%M:%S")
  date2<-strptime(paste(R2$UTC_Date[i+1],R2$UTC_Time[i+1],sep=" "),format="%m/%d/%Y %H:%M:%S")
  dt<-difftime(date2,date1,units="secs")
  R2$deltat_seconds[i+1]<-as.numeric(dt)
  R2$speed<-R2$dist2/as.numeric(dt)
}
proc.time()-start
saveRDS(R2,file="R2.data.rds") 
##############################################################################################################
####thresholding 
#############################################################################################################

thresh<-2000
subsetbelowthreshold<-R2[R2$dist2<thresh,]
length(R2[,1])-length(subsetbelowthreshold[,1])
length(subsetbelowthreshold[,1]) 
plot(subsetbelowthreshold$Longitude, subsetbelowthreshold$Latitude)# misplaced relocation points still in the data even after thresholding  

saveRDS(subsetbelowthreshold,file="R2.data.rds.cleaned") 
Round2 <- readRDS("R2.cleaned.rds")

###################################################################################################################
##################################################################################################################
### create an object of class trajectory for round 2 data
##################################################################################################################
#################################################################################################################
#entire data set for round 2 is saved in an object called R2
locs <- R2$place
locs <-as.data.frame(getutm(R2$Longitude,R2$Latitude)) 
head(locs)
tail(locs)
#points(d1$place)

head(locs)
plot(locs)
plot(locs, pch=20, col=R2$groupcol)

#######First the date needs to be transformed into an object of the class POSIXct. 

class(R2$UTC_Date)#factor
da<-strptime(paste(R2$UTC_Date,R2$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S")  
da <- as.POSIXct(da, tz="UTC") 

head(da)
tail(da)
class(da)

# create an object of class ltraj to store herd movements data: 

tr2 <- as.ltraj(coordinates(locs), date = da, id = R2$idcol)#non unique dates for a given burst.

#find out duplicates
unique(R2$idcol)
for (id in unique(R2$idcol)){
  
  dsubset<-da[R2$idcol==id]
  print(sprintf("%s : %i dates/%i unique",id,length(dsubset),length(unique(dsubset))))
}
summary(duplicated(da))
testdf<-data.frame(da,duplicated(da),R2$idcol)
testdf[testdf$duplicated.da.==T,]
plot(da)

#explore the trajectory data 
class(tr1)
head(tr1[[1]])
tail(tr1[[1]])

head(tr1[[2]]) 
##graphical display of the bursts
plot(tr1)#returns an error
tr.round1 <- ld(tr1) #converting trajectory to data frame
head(tr.round1) #an excellent dataframe 

saveRDS(tr.round1,file="tr.round1.rds")  
dl(tr.round1)#converting a data frame to a trajectory type object. 

is.regular(tr1)#false

#write summary to excel file 
summary(tr1)
install.packages("xlsx")
library(xlsx)
write.xlsx(a, "trajectory_round1_summary.xlsx")#write a trajectory summary in a spread sheet. 
####################################################################################################################
###################################################################################################################
###############Cutting a burst into several segments 
##################################################################################################################
##################################################################################################################




#################################################################################################################
##############################################################################################################
### working with all rounds 
################################################################################################################
#################################################################################################################


m1<-readRDS("cattleGPSfiles")#to read a .RDS file
m2 <- m1[, c(1:3, 5:12, 22,23)]#subset

head(m2)
m2$Longitude<-as.numeric(m2$Longitude)
m2$Elevation<-as.numeric(m2$Elevation)


sum(is.na(m2$Longitude))
sum(is.na(m2$Elevation))
sum(is.na(m2$Latitude))

m2 <- na.omit(m2)
#m2 <- na.omit(m2, cols=c("Longitude", "Latitude", "Elevation", "Speed", "UTC_Date", "UTC_Time", "groupcol", "idcol"))## o
length(m2$Longitude)
length(m2$Latitude) 
length(m2$Elevation)
length(m2$UTC_Date)
length(m2$UTC_Time)
length(m2$idcol)

#new_DF <- k1[is.na(k1$Longitude),]
#head(new_DF)

is.factor(m2$groupcol)
is.factor(m2$idcol)


###################################################################################################################
##################################################################################################################
### create an object of class trajectory
##################################################################################################################
#################################################################################################################
#d1
head(m2)
m2$place <- getutm(m2$Longitude,m2$Latitude)
#points(d1$place)
plot(m2$place)

###############################################################################################################
##############################################################################################################
#######First the date needs to be transformed into an object of the class POSIXct. 
#############################################################################################################
#############################################################################################################
class(m2$UTC_Date)#factor

da<-strptime(paste(m2$UTC_Date,m2$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") #date conversion by mario 
da <- as.POSIXct(da, tz="UTC") 


#da <- as.character(d1$UTC_Date)
head(da)
tail(da)
class(da)

locs <- m2$place
head(locs)
#We can then create an object of class ltraj to store the herd movements:
tr.m2 <- as.ltraj(coordinates(locs), date = da, id=m2$idcol)#Error in as.ltraj(coordinates(locs), date = da, id = m2$idcol) : non unique dates for a given burst


##################duplicate debugging
setwd("D:/shibia/herds/")
data<-readRDS("Round2.rds")
da<-strptime(paste(data$UTC_Date,data$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") 
dapos <- as.POSIXct(da, tz="UTC")
da[1]
dapos[1]
summary(duplicated(da))
summary(duplicated(dapos))
locs <-as.data.frame(getutm(data$Longitude,data$Latitude)) 

##checking for duplicates 
a1<-data[data$idcol=="E038",]
da1<-da[data$idcol=="E038"]
locs1<-locs[data$idcol=="E038",]
length(da1)
length(unique(da1))
summary(duplicated(da1))
dups<-duplicated(da1)
plot(locs1$lon[!dups],locs1$lat[!dups],col="green")
points(locs1$lon[dups],locs1$lat[dups],col="red")
a1[dups,]
?duplicated

dups[dups==T]
which(dups==T)
a1[19890:19899,]
a1[21500:21505,]
a1[22750:22751,]

#remove all duplicates including dates duplicated for different animals-not a good thing though! 
dups<-duplicated(da)
summary(dups)
datawithoutdupes<-data[-which(dups),]
locswithoutdupes<-locs[-which(dups),]
dawithoutdupes<-da[-which(dups)]

#save data in an object of as.ltraj -in this data, even the dates duplicated for differerent animals are dropped
require(adehabitatLT)
tr.m2 <- as.ltraj(coordinates(locswithoutdupes), date = dawithoutdupes, id=datawithoutdupes$idcol)#
length(locswithoutdupes$lon)
length(dawithoutdupes)
length(datawithoutdupes$idcol)
length(unique(dawithoutdupes))


summary(duplicated(dawithoutdupes))
row.names((locswithoutdupes))

######################################################################################################
for (id in unique(data$idcol)){
  dsubs<-da[data$idcol==id]
  print(sprintf("%s: %i/%i: %i dupes!",id,length(dsubs),length(unique(dsubs)),sum(duplicated(dsubs))))
}
#######################################################################################################

dupes_total<-0
data_wodupes<-data[1,]#define a dataframe to later fill in with values 
for (id in unique(data$idcol)){
  dsubs<-da[data$idcol==id]
  datasubset<-data[data$idcol==id,]
  dupes<-which(duplicated(dsubs))
  print(sprintf("%s: %i/%i: %i dupes!",id,length(dsubs),length(unique(dsubs)),length(dupes)))
  dupes_total<-dupes_total+sum(length(dupes))
  data_wodupes<-rbind(data_wodupes,datasubset[-dupes,])
}
print(sprintf("%i total dupes",dupes_total))
dim(data_wodupes)
length(data$Longitude)-length(data_wodupes$Longitude)
length(data_wodupes$Longitude)+416
##################################
#an alternative way yet useful in droping duplicates. date should uniquely define an animal and not duplicated!  

animal_date<-paste(as.character(data$idcol),as.character(da))
summary(duplicated(animal_date))#1842416 false and 416 true
data_wodupes<-data[-which(duplicated(animal_date)),]#less duplicated observation and left the non duplicated data
#check for duplicates again. 
length(data_wodupes$Longitude)+416
length(data$Longitude)-length(data_wodupes$Longitude)


dupes<-which(duplicated(animal_date))#index of dupes
data_wodupes<-data[-dupes,]
locswithoutdupes<-locs[-dupes,]
dawithoutdupes<-dapos[-dupes]
require(adehabitatLT)
tr.m2 <- as.ltraj(coordinates(locswithoutdupes), date = dawithoutdupes, id=data_wodupes$idcol)#

#check out on the lengths -should be of the same the same length
length(locswithoutdupes$lon)
length(dawithoutdupes)
length(data_wodupes$idcol)
#length(unique(data_wodupes))

#explore the trajectory data 
class(tr.m2)
head(tr.m2[[1]])
tail(tr.m2[[1]])

head(tr.m2[[2]]) 
##graphical display of the bursts
#plot(tr.m2)#returns an error
tr.round2 <- ld(tr.m2) #converting trajectory to data frame
head(tr.round2) #an excellent dataframe 

saveRDS(tr.round2,file="tr.round2.df.rds")  
dl(tr.round2)#converting a data frame to a trajectory type object. 

is.regular(tr.m2)#false

# #write summary to excel file 
# summary(tr1)
# install.packages("xlsx")
# library(xlsx)
# write.xlsx(a, "trajectory_round1_summary.xlsx")#write a trajectory summary in a spread sheet.


#hunting of the duplicates in round 2 once again. 
id<-"E038"
ids<-which(data$idcol==id)#now index ids in the data 
summary(duplicated(da[ids]))
which(duplicated(da[ids])) #in the data

range(ids) #what is the range of the index
which(duplicated(data[data$idcol==id,]))#know IDs containing duplicated date or data

id<-"E028"
ids<-which(data$idcol==id)
range(ids)
summary(duplicated(da[ids]))
which(duplicated(da[ids]))

##################################################################################################################
#################################################################################################################
######## round 3
#################################################################################################################
#################################################################################################################

##setwd to data folder 
list.files()
list_files <- list.files()
length <- length(list_files)

##write cattleTracks to R data frame data type
###identiy a group and animal
start <- proc.time()
newdata<-data.frame()# pre-allocate some variable space for our output data
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

#setwd
head(newdata)
R3 <- newdata
saveRDS(newdata,file="Round3.rds") 

n.data <- readRDS("Round3.rds")
n.data.2 <- n.data[n.data$Longitude=="NA",]
head(n.data.2)
#new_DF <- DF[is.na(DF$Var),] #an example

d2 <- newdata
sum(is.na(newdata$Longitude))
sum(is.na(newdata$Latitude))
sum(is.na(newdata$Elevation))

sum(is.na(newdata$UTC_Date))
sum(is.na(newdata$UTC_Time))

# newdata$Longitude<-as.numeric(newdata$Longitude)
# newdata$Latitude <- as.numeric(newdata$Latitude)
# newdata$Elevation<-as.numeric(newdata$Elevation)


d2<- d2[, c(1:12, 22,23)]#subset 
d2<- na.omit(d2)##create a new data frame without a missing value. 

sum(is.na(d2$Longitude))
sum(is.na(d2$Latitude))
sum(is.na(d2$Elevation))


length(newdata$Longitude)-length(d2$Longitude)

da<-strptime(paste(d2$UTC_Date,d2$UTC_Time,sep=" ", tz="UTC"),format="%m/%d/%Y %H:%M:%S") 
dapos <- as.POSIXct(da, tz="UTC")
da[1]
dapos[1]

summary(duplicated(da)) #includes date duplicated for different animals.  
summary(duplicated(dapos))#includes date duplicated for different animals 

locs <-as.data.frame(getutm(d2$Longitude,d2$Latitude)) 
head(locs)
tail(locs)
sum(is.na(locs$lon))#na present 0

### only drop duplicate observations for a single ID excluding duplication for different animals.
## dig deeper in searching for duplicates; http://www.cookbook-r.com/Manipulating_data/Finding_and_removing_duplicate_records/  

# Show the repeat entries
#data_widupes <- df[duplicated(df),] #find data with dupes, an example 

animal_date<-paste(as.character(d2$idcol),as.character(da)) #paste together two characters of different length. 
summary(duplicated(animal_date))## Is each row a repeat? (returns summary stats 1157412 false and 0 true duplicates)  
data_wodupes<-d2[-which(duplicated(animal_date)),]#data without duplicates**** returns 0, why?
data_wodupes <- d2

#confirmation of whether the hunt for duplicates is correctly executed. 
length(data_wodupes$Longitude)+0 #should give us the correct length of longitude in the original dataset d2. 
length(d2$Longitude)-length(data_wodupes$Longitude)#give us the number of duplicated entries 


dupes<-which(duplicated(animal_date))#index of dupes
data_wodupes<-d2#[-dupes,]
locswithoutdupes<-locs#[-dupes,]
dawithoutdupes<-dapos#[-dupes]
require(adehabitatLT)
tr.m3 <- as.ltraj(coordinates(locswithoutdupes), date = dawithoutdupes, id=data_wodupes$idcol)#

#check out on the lengths -should be of the same the same length
length(locswithoutdupes$lon)
length(dawithoutdupes)
length(data_wodupes$idcol)
#length(unique(data_wodupes))

#explore the trajectory data 
class(tr.m3)
head(tr.m3[[1]])
tail(tr.m3[[1]])

head(tr.m3[[2]]) 
##graphical display of the bursts
#plot(tr.m2)#returns an error
tr.round3.df <- ld(tr.m3) #converting trajectory to data frame
head(tr.round2) #an excellent dataframe 

saveRDS(tr.round3.df,file="tr.round3.df.rds")  
#dl(tr.round3.df)#converting a data frame to a trajectory type object. 
head(tr.round3.df)
tr.round3.df$speed<-tr.round3.df$dist/tr.round3.df$dt

## thresholding dist
thresh<-2000
subsetbelowthreshold<-tr.round3.df[tr.round3.df$dist<thresh,]
length(tr.round3.df[,1])-length(subsetbelowthreshold[,1])
length(subsetbelowthreshold[,1]) 
plot(subsetbelowthreshold$x, subsetbelowthreshold$y)

plot(subsetbelowthreshold$dist, subsetbelowthreshold$speed)
boxplot(subsetbelowthreshold$dist) 

hist(subsetbelowthreshold$dist, breaks=20, freq=FALSE, xlab="step length", main="Histogram of cattle step lengths")

###regularize a trajectory by placing missing values 

head(subsetbelowthreshold$date) 
min(subsetbelowthreshold$date, na.rm = TRUE)
max(subsetbelowthreshold$date, na.rm = TRUE)

#plotltr(tr.m3[[1]], "dt/3600/24") 

refda <- strptime("2013-09-02 00:00:00", "%Y-%m-%d %H:%M:%S") #define a reference date and replace missing values with na  
tr.m3.test <- setNA(tr.m3, refda, 5, units = "min")
head(tr.m3.test[[1]])
is.regular(tr.m3.test)#FALSE

refda <- strptime("00:00:00", "%H:%M:%S", tz="UTC") #now set time step to the actual time step
tr.m3.test.2 <- sett0(tr.m3.test, refda, 5, units = "min")

is.regular(tr.m3.test.2)#TRUE
############################################################################################################
###### set limits 
############################################################################################################

is.sd(tr.m3.test.2) #is it trajectory of same duration  and returns FALSE

##We can use the function set.limits to define the time of beginning and ending
#of the trajectories. This function adds NAs to the beginning and ending of the monitoring when required: 

# tr.m3.test.2 <- set.limits(ib3, begin = "2003-06-01 00:00",
#                   dur = 14, units = "day", tz="Europe/Paris", pattern = "%Y-%m-%d %H:%M") #set limit to same duration

rd.tr.m3.test.2<- rdSteps(tr.m3.test.2, reproducible = TRUE) #generate 10 random steps for every presence data 
head(rd.tr.m3.test.2)

saveRDS(rd.tr.m3.test.2, file = "rd.tr.m3.test.2.df.rds")# a data frame for round 3 data

boxplot(rd.tr.m3.test.2$rel.angle~rd.tr.m3.test.2$case)
boxplot(rd.tr.m3.test.2$dist~rd.tr.m3.test.2$case)




