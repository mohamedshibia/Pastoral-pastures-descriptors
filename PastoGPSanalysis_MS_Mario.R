
#setwd
library(sf)
#setwd("~/mnt/agroscope_os/2/2/5/2/3/1935/Shibia")
setwd("H:/Fern153/Gompertz/PastoGPS/30082018")
#Load functions
source("PastoGPS functions.R")

#Retrieve stored data
res_all_df <- readRDS("res_all.rds")


##analyse movement parameters
###a bit of data cleaning and analysis 

sel <- which(res_all_df$asym<20000 & res_all_df$asym>0 & !is.na(res_all_df$asym))
res_all_df2 <- droplevels(res_all_df[sel,]) ###drop outliers 

##Calculating median values for individuals and month within year 
res_all_df2$month <-  format(as.Date(as.character(res_all_df2$date), format="%d.%m.%y"),"%m")
res_all_df2$year <-  format(as.Date(as.character(res_all_df2$date), format="%d.%m.%y"),"%Y")
res_all_df2$dayid <- format(as.Date(as.character(res_all_df2$date), format="%d.%m.%y"),"%j_%Y")#doy

#table(res_all_df2$dayid,res_all_df2$id)

#begin1
#boxplots for all individuals  
pdf("Figure asym vs id + month.pdf")
boxplot(asym ~ month+id, data=res_all_df2)
dev.off()

#output in begin1 is not elegant although we could still find seasonal differences in asym (asymptote) among id (animal id). I would think of how to improve output in begin1. 


#one way is to separate boxplots by id so that we have 22 of them (option1)
#most complicated but intuitive way is to look at the data on the basis of dayid -a variable for day of year (option2) - 
#or any other interesting way you find suitable

png("H:/Fern153/Gompertz/PastoGPS/19072018/test.png",height=120,width=30,units="cm",res=300)
par(mfrow=c(11,2))
for (id_tmp in unique(res_all_df2$id)){
  data_tmp <- res_all_df2[res_all_df2$id==id_tmp,]
  boxplot(asym ~ month, data=data_tmp,main=id_tmp,xlab="Months", ylab="Asymptote",ylim=c(0,20000))
}
dev.off()

##An alternative presentation of grazing lenth-suggested by MARIO
png("F:/Fern153/Gompertz/PastoGPS/30082018/test_gl.png",height=120,width=30,units="cm",res=300)
par(mfrow=c(11,2))
for (id_tmp in unique(res_all_df2$id)){
  data_tmp <- res_all_df2[res_all_df2$id==id_tmp,]
  boxplot(GrazLength ~ month, data=data_tmp,main=id_tmp,xlab="Months", ylab="Grazing length",ylim=c(0,1400))
}
dev.off()

#lineplot of monthly medians
#aggregate(res_all_df2$GrazLength,by=list("year","month"),FUN = median)

grazlen_monthly <- aggregate(GrazLength~year+month,data=res_all_df2,FUN = median)
grazlen_monthly$month <- as.numeric(grazlen_monthly$month)
plot(c(min(grazlen_monthly$month),max(grazlen_monthly$month)),c(0,max(grazlen_monthly$GrazLength)),type="n")
for (year in unique(grazlen_monthly$year) ){
  subs <-grazlen_monthly[grazlen_monthly$year==year,]
  lines(subs$month,subs$GrazLength)
  text(subs[1,2],subs[1,3],label=subs[1,1])
}

#lineplots of means 

grazlen_mean <- aggregate(GrazLength~year+month,data=res_all_df2,FUN = mean)
grazlen_mean$month <- as.numeric(grazlen_mean$month)
plot(c(min(grazlen_mean$month),max(grazlen_mean$month)),c(0,max(grazlen_mean$GrazLength)),type="n")
for (year in unique(grazlen_mean$year) ){
  subs <-grazlen_mean[grazlen_mean$year==year,]
  lines(subs$month,subs$GrazLength)
  text(subs[1,2],subs[1,3],label=subs[1,1])
}

#lineplots of asym


grazlen_asym <- aggregate(asym~year+month,data=res_all_df2,FUN = median)
grazlen_asym$month <- as.numeric(grazlen_asym$month)
plot(c(min(grazlen_asym$month),max(grazlen_asym$month)),c(0,max(grazlen_asym$asym)),type="n")
for (year in unique(grazlen_asym$year) ){
  subs <-grazlen_asym[grazlen_mean$year==year,]
  lines(subs$month,subs$asym)
  text(subs[1,2],subs[1,3],label=subs[1,1])
}


# res_all_df2$yearmonth <- paste(res_all_df2$month, res_all_df2$year)
# tapply(res_all_df2$asym, list(res_all_df2$id, res_all_df2$yearmonth), median)

setwd("H:/Fern153/Gompertz/PastoGPS/30082018")
res_all_df3 <- readRDS("res_all_df3.rds") ##contains range units walked at peak grazing hours and still containing extreme asymptotes ********************

#### establish ties and row max on proportion of rru walked during peak hours 



colnums <- (ncol(res_all_df3)-6) : ncol(res_all_df3)
n_ties <- vector()
col_max <- vector()
i <- 1
for (i in 1:nrow(res_all_df3)){
  subs <- res_all_df3[i,]
  subs[,colnums]
  n_ties[i] <- sum(subs[,colnums]==max(subs[,colnums]))
  if (n_ties[i] == 1){
    col_max[i] <- colnames(subs[,colnums])[subs[,colnums]==max(subs[,colnums])]
  }
  else {
    col_max[i] <- NA
   }
}
test <- cbind(res_all_df3[,colnums],n_ties,col_max)
summary(as.factor(test$n_ties))

saveRDS(test, "propRRUmax.rds")#max percent share 
test <- readRDS("propRRUmax.rds")
summary(as.factor(test$n_ties))

# sel2 <- which(test$n_ties>6)#six or more ties  
# test2 <- droplevels(test[sel2,])
# table(test2$n_ties)
# 
# sel3 <- which(!test$n_ties>6)#six or more ties  
# test3 <- droplevels(test[sel3,])
# table(test3$n_ties)

## 
res_all_df3 <- cbind(res_all_df3,n_ties,col_max)
res_all_df3$month <-  format(as.Date(as.character(res_all_df3$date), format="%d.%m.%y"),"%m")
res_all_df3$year <-  format(as.Date(as.character(res_all_df3$date), format="%d.%m.%y"),"%Y")
res_all_df3$dayid <- format(as.Date(as.character(res_all_df3$date), format="%d.%m.%y"),"%j_%Y")#doy
saveRDS(res_all_df3, "res_all_df4.rds") ##number of ties, col max  **************************

#a bit of data cleaning 

dt <- readRDS("res_all_df4.rds") #contains # of ties, max rru in each daily trip  
#nothunt<-dt[dt$n_ties=="1",] #select no ties observations 
max(dt$asym, na.rm = TRUE)
min(dt$asym, na.rm = TRUE)

sel <- which(dt$asym<20000 & dt$asym>0 & !is.na(dt$asym)) #data cleaning 
#res_all_df6 <- droplevels(test3[sel,])
dt2 <- droplevels(dt[sel,]) # extreme observations dropped 

saveRDS(dt2, "res_all_df5.rds")##exclude all extreme observations, and contains ties >6  

################################################
##k-means clustering  daily trips into respective stress levels 
###############################################
#12.12.2018

setwd("H:/Fern153/Gompertz/PastoGPS/30082018") 
res_all_df5 <- readRDS("res_all_df5.rds") 
res_all_df5_ <- na.omit(res_all_df5[,c(1:2, 12:16, 22:23)]) 


ref <- na.omit(res_all_df5_[,c(3:9)]) #less length than original df  

##kmeans 

wss <- (nrow(ref))*sum(apply(ref,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(ref,iter.max = 300, algorithm="MacQueen",
                                     centers = i)$withinss)

plot(1:15, wss, type="b", xlab="Number of clusters", ylab="Within group sum of squares")

set.seed(1234)
km <-kmeans(ref,5, iter.max = 200) 
if (km$ifault==4) { km = kmeans(ref,iter.max = 200, km$centers, algorithm="MacQueen"); }

attributes(km)
km$centers
km$size  
#km$cluster

#aggregate(ref, by=list(km$cluster), FUN=mean)
#plot(bouts.sub$dist,bouts.sub$rel.angle,col=c("red","green","blue","yellow")[km$cluster])
#summary(km$cluster) 

      # asym       b1       b2        k1        k2 TripLength GrazLength
# 1  2439.498 301.6586 821.8257  43.07415  74.44944   799.0127  329.32295
# 2 15037.358 445.3052 657.6431 159.20565 127.57978   931.5321   47.81193
# 3  8446.356 350.6745 759.4108  84.87160  97.83697   940.1371  107.11483
# 4  5110.545 325.3066 785.5761  64.43245  84.68865   867.3044  186.78979

#append cluster assignment 
length(res_all_df5_$asym)  #data frame length              #
length(km$cluster)                                #
res_all_df5_2<- data.frame(res_all_df5_, km$cluster) #

saveRDS(res_all_df5_2, "res_all_df5_2.rds") ##contains clusters 

#res_all_df5.rds is the main df 
#res_all_df5_2 is subset of res_all_df5

install.packages("data.table") 
library(data.table)

a <- data.table(res_all_df5) #original df 
setkey(a, "burst")

b <- data.table(res_all_df5_2) #newly created df 
setkey(b, "burst")

c <- a[b] #select subset within original df 
table(c$km.cluster)

#subsetting rows using i 
c <- readRDS("00kmeansDFdata.rds") ##possible to load data from a disk 
mean(c[c$km.cluster==4,]$asym)
c <- data.table(c)

e <- c[km.cluster==4]
f <- c[km.cluster==3] #waterpoints
g <- c[km.cluster==2] #water points
h <- c[km.cluster==1]

saveRDS(e, "01kmclass4.rds")#3449
saveRDS(f, "01kmclass3.rds")#1794
saveRDS(g, "01kmclass2.rds")#436
saveRDS(h, "01kmclass1.rds")#4493

k <- c[km.cluster %in% c(1,2,3)] #select all rows that have the value 1,2,3 in column km.cluster 

#plotting using pastoGPS functions 

d <- data.frame(c)
saveRDS(d, "00kmeansDFdata.rds")

c <- readRDS("00kmeansDFdata.rds")
#saveRDS(d, "res_all_df7.rds") ## clustering  
##end kmeans 
#require(pastorGPS)
table(c$km.cluster)
pdf("1gr2.pdf", paper="a4r")
for(i in which(c$km.cluster==2)) {
  plotTrip(c[i,])
}
dev.off()

##task plot on a pdf document plots with no rru intercepted 

res_all_df5 <- readRDS("res_all_df5.rds")
table(res_all_df5$n_ties)

pdf("1delNon-rru.pdf", paper="a4r")
for(i in which(res_all_df5$n_ties==7)) {
  plotTrip(res_all_df5[i,])
}
dev.off()

##drop n_ties==7
sel3 <- which(!res_all_df5$n_ties>6)#six or more ties  
test3 <- droplevels(res_all_df5[sel3,])
table(test3$n_ties)
saveRDS(test3, "res_all_df6.rds")# drop observations with six and more ties 

res_all_df6 <- readRDS("res_all_df6.rds") ##no bursts with 0 intercepted rrus

##looking for a suitable bins to create waterpoints buffers, and looks like 0-20km in a step of 5
max(res_all_df6$asym, na.rm = TRUE)
min(res_all_df6$asym, na.rm = TRUE)
fivenum(res_all_df6$asym)# 60.24492  2588.27340  4086.06018  6359.51470 19998.46501

max(res_all_df6$maxDis)
max(test3$maxDis, na.rm = TRUE)
min(test3$maxDis, na.rm = TRUE)
fivenum(test3$maxDis)

range(res_all_df6$maxDis)
hist(res_all_df6$maxDis) 
hist(res_all_df6$asym)
##buffering water points use bins 0-5km; 0-10km; 0-15km; 0-20km

##lm fitting 
##check if maxdis and 

c018.lm = lm(res_all_df6$asym ~ res_all_df6$maxDis, data=res_all_df6) 

#summary(c001.lm) 
summary(c018.lm)$r.squared  
summary(c018.lm)$adj.r.squared
coef(summary(c018.lm))[2,4] 
coef(summary(c018.lm))[2,1] #slope 

coeff=coefficients(c018.lm) 


#par(mgp=c(2,1,0), mar=c(3,3,1,1)) 
eq = paste0("y = ", round(coeff[2],1), "*x ", round(coeff[1],1)) # equation of the line : 
plot(res_all_df6$asym, res_all_df6$maxDis, main=eq) # a nice plot
abline(c018.lm, col="red")


# Scatterplot ; test for pair correlation: rescale data first! 

pairs(dt[c("asym","b1","b2","k1","k2","GrazLength","TripLength")], main="My Title ", pch=22,
      bg=c("red", "yellow", "blue")[unclass(dt$month)]) 

require(caret)
predictors <- dt[, c("asym","b1","b2","k1","k2","GrazLength","TripLength")]

trans = preProcess(predictors, 
                   c("BoxCox", "center", "scale"))
predictorsTrans = data.frame(
  trans = predict(trans, predictors))

pairs(predictorsTrans, main="My Title ", pch=22,
      bg=c("red", "yellow", "blue")[unclass(dt$month)]) 


