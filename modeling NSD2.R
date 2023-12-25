

####script name; modeling NSD2

library(sf)
setwd("E:/Fern153/Gompertz")
#setwd("C:/HomeOffice/201609/Shibia")
datansd1 <- d<-  droplevels(readRDS("data.subset.rds"))


pdf("bursttest.pdf")
for(i in 1:26){
layout(matrix(c(1,1,2,3), 2,2,byrow=F), width=c(2,1))
sub <- d[as.numeric(d$burst)==i,]
sub <- sub[!is.na(sub$x),]
pts <- st_as_sf(sub, coords = c("x", "y"))
lns <- st_linestring(as.matrix(sub[c("x", "y")]))
m_y<-mean(st_bbox(lns)[c(2,4)])

sub$dis <- as.numeric(st_distance(pts[1,], pts))^2 #Calculate distance from 1st point
sub$minsday <- as.POSIXlt(sub$date)$hour*60+as.POSIXlt(sub$date)$min #Calculate mins per day

plot(lns, ylim=m_y+c(3000,-3000), axes=T)
title(paste(sub$burst[1],":", sub$Trip[1], ":", format(sub$date[1], format="%d.%m.%y")))

#Double exponential
M<-try(nls(dis ~  asym /(1+exp((xmidA-minsday)/scale1)) + 
           (-asym /(1 + exp((xmidB-minsday)/scale2))), data=sub, start = c(asym=5000^2,xmidA=300,xmidB=800,scale1=50,scale2=70)), silent=T) #Try is used to continue in case of non-convergence
#Double exponential with upper bound
if(class(M) == "try-error") M1<-try(nls(dis ~  asym /(1+exp((xmidA-minsday)/scale1)) + 
           (-asym /(1 + exp((xmidB-minsday)/scale2))), data=sub, start = c(asym=5000^2,xmidA=300,xmidB=800,scale1=50,scale2=70), upper=c(asym=max(sub$dis)), algorithm="port"), silent=T)
#Double Gompertz
M2<-try(nls(dis ~   asym * (exp(-exp((b1-minsday)/k1)) - exp(-exp((b2-minsday)/k2))), data=sub, start = c(asym=5000^2,b1=300,b2=800,k1=80,k2=5), control=list(maxiter=5000,minFactor=1/1000000)), silent=T)

if(class(M) != "try-error") legend("topleft", legend=round(coef(M),0))
if(class(M) == "try-error" & class(M1) != "try-error") legend("topleft", legend=round(coef(M1),0))
if(class(M2) != "try-error") legend("topright", legend=round(coef(M2),0))

plot(sub$date, sub$dis , t="l",ylim=c(0,6000^2), ylab="Distance from start")
if(class(M) != "try-error") lines(sub$date, fitted(M), col="red", lwd=2)
if(class(M) == "try-error" & class(M1) != "try-error") lines(sub$date, fitted(M1), col="pink", lty=2, lwd=2)
if(class(M2) != "try-error") lines(sub$date, fitted(M2), col="blue")

plot(sub$date, c(0,cumsum(st_distance(pts$geometry[-nrow(pts)], pts$geometry[-1], by_element=T))), t="l", ylim=c(0,20000), ylab="Commulated distance")
}
dev.off()

# Entire dataset (Reads rds for each individual but can easily be adapted for the complete rds)
d <- readRDS("Traj_c001.rds")
ind <- levels(d$id)

out <- do.call("rbind", lapply(ind, function(f){
  d <- readRDS(paste0("Traj_",f,".rds"))
  d <- d[!is.na(d$x),] #Omit records without positions
  res <- do.call("rbind", lapply(levels(droplevels(d$burst)), function(i){
    sub <- d[d$burst==i,]
    sub <- sub[!is.na(sub$x),]
    pts <- st_as_sf(sub, coords = c("x", "y"))
    sub$dis <- as.numeric(st_distance(pts[1,], pts)) #Calculate distance from 1st point
    sub$minsday <- as.POSIXlt(sub$date)$hour*60+as.POSIXlt(sub$date)$min #Calculate mins per day
    #Double exponential
    M<-try(nls(dis ~  asym /(1+exp((xmidA-minsday)/scale1)) + 
           (-asym /(1 + exp((xmidB-minsday)/scale2))), data=sub, start = c(asym=5000,xmidA=300,xmidB=800,scale1=50,scale2=70)), silent=T) #Try is used to continue in case of non-convergence
    Msq<-try(nls(dis^2 ~  asym /(1+exp((xmidA-minsday)/scale1)) + 
           (-asym /(1 + exp((xmidB-minsday)/scale2))), data=sub, start = c(asym=5000^2,xmidA=300,xmidB=800,scale1=50,scale2=70)), silent=T) #Try is used to continue in case of non-convergence
    #Double exponential with upper bound
  M1<-try(nls(dis ~  asym /(1+exp((xmidA-minsday)/scale1)) + 
           (-asym /(1 + exp((xmidB-minsday)/scale2))), data=sub, start = c(asym=5000,xmidA=300,xmidB=800,scale1=50,scale2=70), upper=c(asym=max(sub$dis)), algorithm="port"), silent=T)
  M1sq<-try(nls(dis^2 ~  asym /(1+exp((xmidA-minsday)/scale1)) + 
           (-asym /(1 + exp((xmidB-minsday)/scale2))), data=sub, start = c(asym=5000^2,xmidA=300,xmidB=800,scale1=50,scale2=70), upper=c(asym=max(sub$dis)), algorithm="port"), silent=T)
    #Double Gompertz
    M2<-try(nls(dis ~   asym * (exp(-exp((b1-minsday)/k1)) - exp(-exp((b2-minsday)/k2))), data=sub, start = c(asym=5000,b1=300,b2=800,k1=80,k2=5), control=list(maxiter=5000,minFactor=1/1000000)), silent=T)
    M2sq<-try(nls(dis^2 ~   asym * (exp(-exp((b1-minsday)/k1)) - exp(-exp((b2-minsday)/k2))), data=sub, start = c(asym=5000^2,b1=300,b2=800,k1=80,k2=5), control=list(maxiter=5000,minFactor=1/1000000)), silent=T)

  if(class(M) == "try-error") M_co <- rep(NA,5) else M_co <- coef(M)
  names(M_co) <- c("M_asym", "M_xmidA", "M_xmidB", "M_scale1", "M_scale2")
  if(class(Msq) == "try-error") Msq_co <- rep(NA,5) else Msq_co <- coef(Msq)
  names(Msq_co) <- c("Msq_asym", "Msq_xmidA", "Msq_xmidB", "Msq_scale1", "Msq_scale2")
  if(class(M1) == "try-error") M1_co <- rep(NA,5) else M1_co <- coef(M1)
  names(M1_co) <- c("M1_asym", "M1_xmidA", "M1_xmidB", "M1_scale1", "M1_scale2")
  if(class(M1sq) == "try-error") M1sq_co <- rep(NA,5) else M1sq_co <- coef(M1sq)
  names(M1sq_co) <- c("M1sq_asym", "M1sq_xmidA", "M1sq_xmidB", "M1sq_scale1", "M1sq_scale2")
  if(class(M2) == "try-error") M2_co <- rep(NA,5) else M2_co <- coef(M2)
  names(M2_co) <- c("M2_asym", "M2_k1", "M2_k2", "M2_b1", "M2_b2")
  if(class(M2sq) == "try-error") M2sq_co <- rep(NA,5) else M2sq_co <- coef(M2sq)
  names(M2sq_co) <- c("M2sq_asym", "M2sq_k1", "M2sq_k2", "M2sq_b1", "M2sq_b2")

    return(cbind(data.frame(id = sub$id[1], burst=i, date=format(sub$date[1], format="%d.%m.%y"), start=format(sub$date[1], format="%H:%M"), end=format(sub$date[nrow(sub)], format="%H:%M")), data.frame(t(M_co)), data.frame(t(Msq_co)), data.frame(t(M1_co)), data.frame(t(M1sq_co)), data.frame(t(M2_co)), data.frame(t(M2sq_co))))
}))
return(res)  
}))
write.table(out, "Parameters.txt")

#A bit of exploratory plotting
d3 <- read.table("Parameters.txt")
d3$date <- as.POSIXct(as.character(d3$date), format="%d.%m.%y")
plot(M2_asym ~ date, data=d3, ylim=c(0,30000), col=rainbow(22)[as.numeric(out$id)], pch=16)
plot(d3$date, d3$M2_k2-out$M2_k1, ylim=c(0,1000), col=rainbow(22)[as.numeric(out$id)], pch=16,cex=.4)
plot(out$M2_asym, out$M2_k2-out$M2_k1)
boxplot(M2_asym ~ id, d3)


#Evaluation of non-convergence
sub <- d[as.numeric(d$burst)==4,]
sub <- sub[!is.na(sub$x),]
pts <- st_as_sf(sub, coords = c("x", "y"))
sub$dis <- as.numeric(st_distance(pts[1,], pts))
sub$minsday <- as.POSIXlt(sub$date)$hour*60+as.POSIXlt(sub$date)$min
plot(sub$date, sub$dis , t="l",ylim=c(0,6000), ylab="Distance from start")

M1<-nls(dis ~   asym * (1 /(1+exp((xmidA-minsday)/scale1)) - 1 /(1 + exp((xmidB-minsday)/scale2))), data=sub, start = c(asym=1000,xmidA=300,xmidB=800,scale1=5,scale2=7), control=list(maxiter=50,minFactor=1/10000, tol=1e-06), upper=c(asym=max(sub$dis)), algorithm="port", trace=T)
lines(sub$date, fitted(M1), col="red")

M2<-nls(dis^2 ~   asym * (exp(-exp((b1-minsday)/k1)) - exp(-exp((b2-minsday)/k2))), data=sub, start = c(asym=5000^2,k1=80,k2=5,b1=300,b2=800), control=list(maxiter=5000,minFactor=1/1000000), trace=T)
lines(sub$date, fitted(M2), col="blue")




##start modelling NSD using nlme 
xyplot(nsdall~minsday|Trip,data=datansd1)  
nls(nsdall~ SSlogis(minsday, Asym, xmid, scal), datansd1) #a self-starter to returns starting values 

#this will help identify parameter estimates for use with nlme
m1<-nls(nsdall ~  asym /(1+exp((xmidA-minsday)/scale1)) + 
          (-asym / (1 + exp((xmidB-minsday)/scale2))), #this part defines eqn 1
        start = c(asym=9.664e+06,xmidA=3.733e+02,xmidB=9.0e+02,scale1=2.784e+01,scale2=2.984e+01)                  
        #these are the starting values for each parameter of the equation 
        ,data=na.omit(datansd1), trace = TRUE)   #this is the data
summary(m1)        #this will print a summary of the converged model

#to graph nsd against time, use:
xyplot(nsdall~minsday|id,data=datansd1)

#now try and model the data including random effects
#start with no variation in the explanatory variable
m3<-nlme(nsdall ~  asym /(1+exp((xmidA-minsday)/scale1)) + 
           (-asym /(1 + exp((xmidB-minsday)/scale2))), #the equation
         fixed = list(asym+xmidA+xmidB+scale1+scale2~1),  #fixed effects
         random= asym ~ 1|id, #random effects: asym varies between individuals
         start = c(asym=1.481e+07,xmidA=4.187e+02,xmidB=8.167e+02,scale1=5.173e+01,scale2=7.715e+01)                  
         #starting vlaues for the parameters in the equation
         ,data=na.omit(datansd1))      #the data
print(AIC(m3))           #this will print the AIC of the converged model


#you can change the random effect structure
m4<-nlme(nsdall ~  asym /(1+exp((xmidA-minsday)/scale1)) + 
           (-asym /(1 + exp((xmidB-minsday)/scale2))),
         fixed = list(asym+xmidA+xmidB+scale1+scale2~1),
         random= asym ~ 1|id/Trip,        #random effects: asym varies between 
         #individuals, and also between trips within a single individual
         start = c(asym=1.481e+07,xmidA=4.187e+02,xmidB=8.167e+02,scale1=5.173e+01,scale2=7.715e+01)
         ,data=na.omit(datansd1))
print(AIC(m4)) 

summary(m4)
#you can show the fitted values
fitted(m4)
#normal probability plots;Can we assume our sample of nsd comes from a population that is Normally distributed?
qqnorm(m4)

qqnorm(residuals(m4))
qqline(residuals(m4))
summary(m4)
#the residuals
plot(m4)



##########################################################################################################
##########################################################################################################
######################## characterizing trips
#######################################################################################################
#######################################################################################################
setwd("E:/Fern153/Gompertz")
bouts <- readRDS("all.trajectory.daily.trips.DF.rds") #raw data that was shared with manuel

datansd <- readRDS("all.trajectory.datansd.final.rds") #the datansd and raw data were merged for personal use 
bouts <- readRDS("bouts2018")

c006 <-bouts[with(bouts, id=="c006"),] 

###a loop to drop all the trips with observations less than required number 
bouts$nobs<-0

animal<-"c001" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c002" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}
animal<-"c003" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c003.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c004" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c004.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}
animal<-"c005" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c006" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c007" #select an animal in the dataframe to include in the loop for plotting daily displacements 
# animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
# plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c008" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c009" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c010" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c010.f" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c011" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c011.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c012" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c012.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c013" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c013.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c014" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c014.f" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c015" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

saveRDS(bouts, file = "bouts2018_2.rds")

hist(bouts$nobs) 
hist(bouts$nobs[bouts$nobs>0]) 
hist(bouts$nobs[bouts$nobs>0 & bouts$nobs<80])
#drop all the trips with < 80 observations 

#########################################################################################
########################################################################################
### calculate the max daily travelled distances  
##################################################################################
#bouts <- bouts.2 
bouts <- readRDS("bouts2018_3.rds")
bouts$NSD.max <- 0

sum(is.na(bouts$R2n))
bouts2 <- bouts[is.na(bouts$R2n),]
bouts3 <- bouts[!is.na(bouts$R2n),]

bouts <- bouts3 #without missing values and complete 4659140 observations. 

animal<-"c001" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}



animal<-"c002" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c003" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c003.w" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c004" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c004.w" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c005" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c006" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}


animal<-"c007" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c008" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c009" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c010" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c010.f" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c011" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c011.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c012" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c012.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c013" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c013.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c014" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c014.f" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c015" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- max(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.max"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

saveRDS(bouts, file = "bouts2018_4.rds")

# bouts <- readRDS("bouts2018_3.rds")
# bouts$NSD.max <- 0
# 
# sum(is.na(bouts$R2n))
# bouts2 <- bouts[is.na(bouts$R2n),]
# bouts3 <- bouts[!is.na(bouts$R2n),]
# 
# bouts <- bouts3
# range(bouts$NSD.max)
range(bouts$nobs)
hist(bouts$nobs) 
hist(bouts$nobs[bouts$nobs>0]) 
hist(bouts$nobs[bouts$nobs>0 & bouts$nobs<80])

#########################################################################
####last NSD value 
#######################################################################
install.packages("dplyr")
library(dplyr)

table(bouts$id)
bouts$NSD.last<-0

animal<-"c001" 
#animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}



animal<-"c002" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c003" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c
for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c003.w" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c004" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c004.w" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c005" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c006" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c007" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c008" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c009" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c010" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c010.f" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c011" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c011.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c012" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c012.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c013" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c013.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c014" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c014.f" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}


animal<-"c015" #select an animal in the dataframe to include in the loop for plotting daily displacements 
# animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
# plot(animal_day_subs$minsday,animal_day_subs$R2n)

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$R2n)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"NSD.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

saveRDS(bouts, file = "bouts2018_5.rds")
saveRDS(bouts, file = "all.trajectory.daily.trips.DF.v2.rds")

bouts <- readRDS("all.trajectory.daily.trips.DF.v2.rds")
ta <- names(bouts)
write.csv(ta, file = "data labels.csv")

boxplot(bouts$dist~bouts$hdaylitr2, xlab="Hours of day (UTC))", ylab="Step lengths (m)", main="Step lengths vs. Time of day")
#large distances walked in 5 minutes-obvious gpd errors e.g., >50k meters
plot(bouts$minsday,bouts$R2n)
#farther places from the start due to errors in GPS
 
###########################
#thresholding of NSD and step lengths  
########################
##missing values are already dropped 
## find out outliers
plot(bouts$minsday, bouts$R2n) # a wise idea is to cut off R2n at 5e+08 
thresh <- 5e+08 #500000000 #that's at about 14km 
bouts.1 <- bouts[bouts$R2n<thresh,] #trancated dataset using NSD value<5e+08


#checks 
length(bouts$x)-length(bouts.1$x)#dropped 87 observations  
length(bouts.1$x)+ 519 #2129929
plot(bouts.1$minsday, bouts.1$R2n)


#thresholding dist at 2000 meters 
thresh <- 2000 
bouts.2 <- bouts.1[bouts.1$dist<thresh,]

#checks 
length(bouts.1$x)-length(bouts.2$x)#dropped 87 observations  
length(bouts.2$x)+ 287 #2129929
plot(bouts.2$minsday, bouts.2$R2n)

boxplot(bouts.2$dist~bouts.2$hdaylitr2, xlab="Hours of day (UTC))", ylab="Step lengths (m)", main="Step lengths vs. Time of day")

saveRDS(bouts.2, file = "all.trajectory.daily.trips.DF.v3.rds")
bouts <- readRDS("all.trajectory.daily.trips.DF.v3.rds")

sum(is.na(bouts$R2n))
bouts2 <- bouts[is.na(bouts$R2n),]#208306
bouts <- bouts[!is.na(bouts$R2n),]#4450028
saveRDS(bouts, file = "all.trajectory.daily.trips.DF.v4.rds")
bouts <- readRDS("all.trajectory.daily.trips.DF.v4.rds") ## data for sharing 

##recode 
library(lattice)
library(car)

#NSD.last recode 
range(bouts$NSD.last)
max(bouts$NSD.last)#38km away
#last(mtcars$mpg)
hist(bouts$NSD.last)

median(bouts$NSD.last)
quantile(bouts$NSD.last)

##a trip ending at the farthest places i.e., when herds are relocated to a different camp. 

# 0%          25%          50%          75%         100% 
# 3.391284e-02 1.176549e+02 5.029228e+02 1.948237e+03 1.518962e+09 

fivenum(bouts$NSD.last)

bouts$trips<- bouts$NSD.last

#################

bouts$NSD_last_rec<- bouts$NSD.last
bouts$NSD_last_rec <- recode(bouts$NSD_last_rec, "3.391284e-02:1.176549e+02='central place'; 1.176549e+02:5.029228e+02='closeby'; 5.029228e+02:1.948237e+03='Fairly away'; 1.948237e+03:1518962499='far away'" )


bouts_nsd.faraway<-bouts[bouts$NSD_last_rec=="far away",] #700-22300 meters away 
plot(bouts_nsd.faraway$minsday, bouts_nsd.faraway$R2n)

hist(bouts_nsd.faraway$NSD.last)
max(bouts_nsd.faraway$NSD.last) # 1518962497 e.g. for distance 38973 meters
range(bouts_nsd.faraway$NSD.last)
min(bouts_nsd.faraway$NSD.last) #1948.237 meters away
median(bouts_nsd.faraway$NSD.last) #5147.148 meters away 


xyplot(R2n~minsday|NSD_last_rec,data=bouts_nsd.faraway)

quantile(bouts_nsd.faraway$NSD.last)# we have problematic trips in upper quantile 

fivenum(bouts_nsd.faraway$NSD.last)

bouts_nsd.faraway$NSD_last_rec2<- bouts_nsd.faraway$NSD.last
bouts_nsd.faraway$NSD_last_rec2 <- recode(bouts_nsd.faraway$NSD_last_rec2, "1.948237e+03:2.875550e+03='central place'; 2.875550e+03:5.147148e+03='closeby'; 5.147148e+03:2.740441e+04='Fairly away'; 2.740441e+04:1518962500='far away'" )


bouts_nsd.faraway$NSD_last_rec2 <- as.factor(bouts_nsd.faraway$NSD_last_rec2)
table(bouts_nsd.faraway$NSD_last_rec2)

boxplot(bouts$NSD_last_rec2, bouts$dist)
xyplot(R2n~minsday|NSD_last_rec2,data=bouts_nsd.faraway)
# 
# ##okay twende kazi tena 
bouts_nsd.faraway2<-bouts_nsd.faraway[bouts_nsd.faraway$NSD_last_rec2=="far away",] 

fivenum(bouts_nsd.faraway2$NSD.last)

bouts_nsd.faraway2$NSD_last_rec3<- bouts_nsd.faraway2$NSD.last
bouts_nsd.faraway2$NSD_last_rec3 <- recode(bouts_nsd.faraway2$NSD_last_rec3, "2.740441e+04:1.270186e+05='central place'; 1.270186e+05:2.060399e+06='closeby'; 2.060399e+06:3.113029e+07='Fairly away'; 3.113029e+07:1518962500='far away'" )


bouts_nsd.faraway2$NSD_last_rec3 <- as.factor(bouts_nsd.faraway2$NSD_last_rec3)
table(bouts_nsd.faraway2$NSD_last_rec3)

boxplot(bouts_nsd.faraway2$NSD_last_rec3, bouts$dist)
xyplot(R2n~minsday|NSD_last_rec3,data=bouts_nsd.faraway2)

###more steps 
bouts_nsd.faraway3<-bouts_nsd.faraway2[bouts_nsd.faraway2$NSD_last_rec3=="far away",] 

fivenum(bouts_nsd.faraway3$NSD.last)

bouts_nsd.faraway3$NSD_last_rec4<- bouts_nsd.faraway3$NSD.last
bouts_nsd.faraway3$NSD_last_rec4 <- recode(bouts_nsd.faraway3$NSD_last_rec4, "31859000:57358245='central place'; 57358245:120827691='closeby'; 120827691:198200298='Fairly away'; 198200298:715489924='far away'" )


bouts_nsd.faraway3$NSD_last_rec4 <- as.factor(bouts_nsd.faraway3$NSD_last_rec4)
table(bouts_nsd.faraway3$NSD_last_rec4)


xyplot(R2n~minsday|NSD_last_rec4,data=bouts_nsd.faraway3)



#########################################################################################
# bouts_nsd.faraway$NSD_last_rec<- bouts_nsd.faraway$NSD.last
# bouts_nsd.faraway$NSD_last_rec <- recode(bouts_nsd.faraway$NSD_last_rec, "30075703:51028691='central place'; 51028691:103556932='closeby'; 103556932:177933520='Fairly away'; 177933520:497328430='far away'" )
# 
# bouts_nsd.faraway3$NSD_last_rec4 <- as.factor(bouts_nsd.faraway3$NSD_last_rec4)
# table(bouts_nsd.faraway3$NSD_last_rec4)
# xyplot(R2n~minsday|NSD_last_rec4,data=bouts_nsd.faraway3)

########################################
###charecterizing trips based on the last NSD. class breaks for trips are informed by categories of last.nsd 

bouts <- readRDS("all.trajectory.daily.trips.DF.v4.rds") ## data for sharing 
bouts$trips <- bouts$NSD.last
bouts$trips <- recode(bouts$trips, "0.0339128373912575:1.176549e+02='return1'; 117.6549:502.9228='return2'; 502.9228:1948.237='return3'; 1948.237:2875.55='return4';2875.55:5147.148='return5';5147.148:27404.41='return6';27404.41:127018.6='return7'; 127018.6:2060399='return8'; 2060399:31130290='return9'; 31130290:57358245='return10';57358245:120827691='return11'; 120827691:198200298='return12'; 198200298:715489924='return13'")

xyplot(R2n~minsday|trips,data=bouts)

bouts$trips <- as.factor(bouts$trips)
table(bouts$trips)

saveRDS(bouts, file = "all.trajectory.daily.trips.DF.v5.rds")


##############################################################
###NSD.max -data grouping 
############################################################

setwd("H:/Fern153/Gompertz")
bouts <- readRDS("all.trajectory.daily.trips.DF.v5.rds")

hist(bouts$NSD.max)
min(bouts$NSD.max)
range(bouts$NSD.max) # 2.143779e+02 2.629026e+09 
median(bouts$NSD.max) #17260743 

###############################################################
#group trips based on nsd max e.g., using the quantiles 

hist(unique(bouts$NSD.max[bouts$NSD.max<500000000]),breaks=50) ##lesss than 22km away
hist(unique(bouts$NSD.max[bouts$NSD.max>500000000]),breaks=50) ##more than 22km away
hist(unique(bouts$NSD.max[bouts$NSD.max>200000000]),breaks=50)##more than 14 km away
hist(unique(bouts$NSD.max[bouts$NSD.max>40000]),breaks=50)##more than 200 m away
hist(unique(bouts$NSD.max[bouts$NSD.max<40000]),breaks=50)##more than 200 m away

dist<-sqrt(unique(bouts$NSD.max))
hist(dist,breaks=1000)
#hist(dist[dist<1000],breaks=100)

windows()
hist(dist[dist<500],breaks=100)
hist(dist[dist>500],breaks=100)

quantile(bouts$NSD.max)
quantile(dist)

##recode 
library(lattice)
library(car)
##farthest point (nsd max) herds ever reached in every day trip  
bouts$place <- bouts$NSD.max
bouts$place <- recode(bouts$place, "2.143779e+02:40000='At home'; 40000:6.627102e+06='nearby'; 6.627102e+06:1.726074e+07='intermediate'; 1.726074e+07:4.050272e+07='far';4.050272e+07:2.629026e+09='very far'")

##drop the following observations
#at home is conotated with the GPS collars staying at home

bouts$place <- as.factor(bouts$place)
table(bouts$place)
boxplot(bouts$place, bouts$dist)
xyplot(R2n~minsday|place,data=bouts)

#miles <- sqrt(bouts$NSD.max)
bouts$nsd2 <- sqrt(bouts$NSD.max)#~max distance (not in nsd) of herd location from the start place

#bouts$burst[which(bouts$nsd2<30)]
bouts$burst[which(sqrt(bouts$NSD.max)<30)] # 

smalldistances_bursts<-bouts$burst[which(sqrt(bouts$NSD.max)<400)]
summary(bouts$id[smalldistances_bursts]) 


# par(mfrow=c(5,5))
# for (i in 1:25){
#   day_animal_subset<-bouts[bouts$burst==smalldistances_bursts[i],]
#   plot(day_animal_subset$x,day_animal_subset$y)
# }
#dev.off()
day_animal_subset<-bouts[bouts$burst=="c001.1278",]
plot(day_animal_subset$x,day_animal_subset$y)

bouts$orbit <-bouts$NSD.max/bouts$NSD.last # ORBIT is defined as a ratio of nsdmax>nsdlast
bouts$orbit_rec <- bouts$orbit>1 # orbit-rec is true for condition meeting the nsdmax>nsdlast.confirm return trips 

table(bouts$orbit_rec)

##subsetting 
longwalks <- bouts[bouts$orbit_rec=="FALSE",] #nsd last and nsd max are the same (5126), trips without a sign of going back to the start place. 

shortwalks <- bouts[bouts$orbit_rec=="TRUE",] #
xyplot(R2n~minsday|orbit_rec, data=longwalks) #drop 
xyplot(R2n~minsday|orbit_rec, data=shortwalks)


boutsub <- shortwalks #there is an evidence of a retun, but not always to the start place.

saveRDS(bouts, file = "all.trajectory.daily.trips.DF.v6.rds")#dataframe complete 
saveRDS(boutsub, file = "all.trajectory.daily.trips.DF.v7.rds")#dropped trips without sign of going back to the start place.i.e., nsdmax=nsdlast 

##task -to filter and drop trips from those either animal or gps collar stayed at home
##thresholding of nsd2
##but first thing first is to identifying cut point

boutsub <- readRDS("all.trajectory.daily.trips.DF.v7.rds")

range(bouts$nsd2)
#range(bouts$NSD.max)

#table(boutsub$place, boutsub$orbit_rec)
table(boutsub$place, boutsub$trips)

##plot outcome to see if the points are stationed at a single place or sign of dispersal

c1 <-boutsub[with(boutsub, place=="nearby" & trips=="return1"),] #subset 
range(c1$nsd2) 
write.csv(c1, file = "nearby_return1.csv")

windows()
hist(c1$nsd2)
hist(unique(c1$nsd2[c1$nsd2<400]),breaks=50) ##lesss than 22km away
dev.off()
##trancate nsd2>800

# c1A2<-boutsub[with(boutsub, place=="very far" & trips=="return13"),]
# range(c1A2$nsd2) 
# write.csv(c1A2, file = "veryfar_return13.csv")
# 
# 
# c1A3<-boutsub[with(boutsub, place=="very far" & trips=="return12"),]
# range(c1A3$nsd2) 
# write.csv(c1A3, file = "veryfar_return12.csv")
# xyplot(R2n~minsday|trips, data=c1A3)
# 
# 
# c1A4<-boutsub[with(boutsub, place=="At home" & trips=="return4"),]
# range(c1A4$nsd2) #51.85089 198.10771-drop
# 
# c1A5<-boutsub[with(boutsub, place=="At home" & trips=="return5"),]
# range(c1A5$nsd2) #58.76102 193.68629-drop
# 
# c1A6<-boutsub[with(boutsub, place=="At home" & trips=="return6"),]
# range(c1A6$nsd2) #58.76102 193.68629-drop
# 
# c1A7<-boutsub[with(boutsub, place=="At home" & trips=="return7"),]
# range(c1A7$nsd2) #58.76102 193.68629-drop
# 
# write.csv(c1A7, file = "athome_return7.csv")
# 
# plot(c1A7$x, c1A7$y, col=c1A7$burst)
# plot(c1A6$x, c1A6$y, col=c1A6$burst)
# plot(c1A5$x, c1A5$y, col=c1A5$burst)
# plot(c1A4$x, c1A4$y, col=c1A4$burst)
# plot(boutsub$x, boutsub$y, col=boutsub$place)
# 
# ####

##farthest point (nsd max) herds ever reached in every day trip  
# boutsub$nsd2recode <- boutsub$nsd2
# boutsub$nsd2recode <- recode(boutsub$nsd2recode, " 14.64165:800='At home'; 800:2574.94532='nearby'; 2574.94532:6364.26109='intermediate'; 6364.26109:51274.02759166='far'")#<800 athome

##         0%         25%         50%         75%        100% 
##14.64165  2574.94532  4157.79586  6364.26109 51274.02759 

# boutsub$nsd2recode <- as.factor(boutsub$nsd2recode)
# table(boutsub$nsd2recode)
# boxplot(boutsub$nsd2recode, boutsub$dist)
# xyplot(R2n~minsday|nsd2recode,data=boutsub)

##subsetting 

thresh <- 800 #thresholding dist at 800 meters 
boutsub.2 <- boutsub[boutsub$nsd2>thresh,] ##no spending at home 
boutsub.3 <- boutsub[boutsub$nsd2<thresh,] ##certainly gps collared herd (or GPS) stayed at home
length(boutsub$x)-length(boutsub.2$x)

saveRDS(boutsub.2, file = "all.trajectory.daily.trips.DF.v8.rds")
saveRDS(boutsub.3, file = "athomegps.rds")

write.csv(boutsub.3, file = "athome.csv")#confirmed at home always
range(boutsub.3$nsd2)

#####################################################################################
####### cleaning GPS collar data
#######################################################################################
setwd("H:/Fern153/Gompertz")
library(lattice)
library(car)

bouts <- readRDS("all.trajectory.daily.trips.DF.v8.rds")


# ##############selecting day events only 
plot(bouts$minsday, bouts$R2n)#herds are quite active as from 200 minutes into the past 1000 minutes of day

range(bouts$minsday)#0 1435 
range(bouts$hdaylitr2)#0 23 
boxplot(bouts$dist~bouts$hdaylitr2, xlab="Hours of day (UTC))", ylab="Step lengths (m)", main="Step lengths vs. Time of day") 
# boxplot(d4$R2n~d4$hdaylitr2.2, xlab="Hours of day (UTC))", ylab="NSD (km2)", main="NSD vs. Time of day")
# 
# ##select daylights hours  
bouts.day <- subset(bouts, bouts$minsday>180 & bouts$minsday<1000) #
plot(bouts.day$minsday, bouts.day$R2n)#herds are quite active as from 200 minutes into the past 1000 minutes of day

bouts.day <- subset(d4, d4$hdaylitr2>3 & d4$hdaylitr2<16) # daylights covers >6Am EAT to < 7.00PM EAT 
saveRDS(bouts.day, file = "all.trajectory.daily.trips.DF.v9.dayevents.rds") #trancated dataset 1000<minsday>200
saveRDS(bouts.day, file = "all.trajectory.daily.trips.DF.v9.b.dayevents.rds") #trancated dataset 1000<minsday>180

boutset <- bouts.day[,1:2]
write.csv(boutset, file = "bouts_dayevents.csv")

ta <- names(bouts.day)
write.csv(ta, file = "data labels_dayevents.csv")


#########################################
############## separate herding trips by camp types 
##################################################

setwd("H:/Fern153/Gompertz") 
library(foreign)
library(lattice)
library(car)
require(rgdal)
require(raster) 
# 
# bouts <- readRDS("all.trajectory.daily.trips.DF.v9.dayevents.rds")
# round1 <- read.csv("DATAROUND1_2011.csv") #supposedly camping sites but returned incorrect points 
# 
# plot(round1$LatDD, round1$lonDD, col=round1$HHID)
# 
# 
# c(round1$lonDD[1],round1$LatDD[1])
# getutm<-function(lon,lat){
#   xy <- cbind(lon, lat)
#   # Convert it to UTM coordinates (in units of meters)
#   return(project(xy, "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))#projected coordinate systems for distance calculations.
# }
# 
# getutm(round1$lonDD[1],round1$LatDD[1]) 
# 
# locs <-as.data.frame(getutm(round1$lonDD,round1$LatDD)) # 
# head(locs)
# plot(locs$lon, locs$lat, col=round1$HHID)
# 
# lon <- locs$lon
# lat <- locs$lat
# 
# round1$lon <- lon
# round1$lat <- lat
# write.csv(round1, file="DATAROUND1_2011_2.csv")
# 
# #shortdist <- bouts$burst[which(sqrt(bouts$NSD.max))<1000]

#########################################
############## separate herding trips by camp types 
##################################################
setwd("H:/Fern153/Gompertz") 
bouts <- readRDS("all.trajectory.daily.trips.DF.v9.dayevents.rds")

round1 <- read.dbf("GPSdayevents_dist.dbf") # wara herds are found within a wara range of 1000 meters. 

library(foreign)
library(lattice)
library(car)
require(rgdal)
require(raster) 
library(dplyr)

bouts$homerange <- round1$RASTERVALU#bouts with the values of 1000 meters are within the range of permanent homesteads  
##data evaluations to identify herds that returned to the start place 
##trips is calculated on the last nsd- indicating extent the trip returned to start place 
##place is calculated on the max nsd -give a clue on the farthest place herd spent from the start place 

### a bit on desribing trips 
table(bouts$place, bouts$trips)

#plotting nsd 

 
animal_day_subs<-bouts[bouts$trips=="return1" & bouts$place=="far",]
plot(animal_day_subs$minsday,animal_day_subs$R2n,type="l", xlab="Minutes of day", ylab="Net squared displacement (NSD)", main="NSD vs Minutes of day") 

animal_day_subs2<-bouts[bouts$trips=="return1" & bouts$place=="intermediate",]
plot(animal_day_subs2$minsday,animal_day_subs2$R2n,type="l", xlab="Minutes of day", ylab="Net squared displacement (NSD)", main="NSD vs Minutes of day") 


animal_day_subs3<-bouts[bouts$trips=="return1" & bouts$place=="nearby",]
plot(animal_day_subs3$minsday,animal_day_subs3$R2n,type="l", xlab="Minutes of day", ylab="Net squared displacement (NSD)", main="NSD vs Minutes of day") 

animal_day_subs4<-bouts[bouts$trips=="return1" & bouts$place=="very far",]
plot(animal_day_subs4$minsday,animal_day_subs4$R2n,type="l", xlab="Minutes of day", ylab="Net squared displacement (NSD)", main="NSD vs Minutes of day") 

xyplot(R2n~minsday|trips, data=bouts)
xyplot(R2n~minsday|place, data=bouts)

# ##################create a new field and assign it a unique id for camp types; 1. base camp and 2. temporary camp
#is the return at the permanent home or temprary camp? 
bouts$camps <- 0
bouts$camps[bouts$trips=="return1" & bouts$homerange==1000] <- 1 # 
bouts$camps[bouts$trips=="return2" & bouts$homerange==1000] <- 1 # 
bouts$camps[bouts$trips=="return3" & bouts$homerange==1000] <- 1 #
bouts$camps[bouts$trips=="return4" & bouts$homerange==1000] <- 1 
bouts$camps[bouts$trips=="return5" & bouts$homerange==1000] <- 1 
bouts$camps[bouts$trips=="return6" & bouts$homerange==1000] <- 1 
bouts$camps[bouts$trips=="return7" & bouts$homerange==1000] <- 1
bouts$camps[bouts$trips=="return8" & bouts$homerange==1000] <- 1 
bouts$camps[bouts$trips=="return9" & bouts$homerange==1000] <- 1
bouts$camps[bouts$trips=="return10" & bouts$homerange==1000] <- 1 
bouts$camps[bouts$trips=="return11" & bouts$homerange==1000] <- 1
bouts$camps[bouts$trips=="return12" & bouts$homerange==1000] <- 1 
bouts$camps[bouts$trips=="return13" & bouts$homerange==1000] <- 1 

table(bouts$camps)
xyplot(R2n~minsday|camps, data=bouts)
# 
saveRDS(bouts, file = "all.trajectory.daily.trips.DF.v10.dayevents.rds")


##summary statistics calculated on daily trips 
library(psych)  
s.bouts <- describeBy(bouts$dist, bouts$burst, mat = TRUE) # one grouping variable
s.bouts2 <-describeBy(bouts$dist, list(bouts$id,bouts$burst), mat = TRUE)#two grouping variables

ta<- table(bouts$dayid,bouts$id)
write.csv(ta, file = "animal_codes_updated20180417.csv")
#table(bouts$id,bouts$id) 



# 
# 
# ##
# 
# 
# #thresholding
# thresh <- 100
# bouts.1 <- bouts[bouts$RASTERVALU>thresh,]
# 
# #checks 
# length(bouts$x)
# length(bouts.1$x)
# length(bouts$x)-length(bouts.1$x)
# 
# bouts1<-bouts[bouts$NSD.max<thresh,]
#long.dry <- bouts.subset[bouts.subset$season=="Long Dry",]#subsetting


##########################################################################################


#########################################################################
####getting the last homerange value -and describing burst on the basis of these parameters
#######################################################################

setwd("H:/Fern153/Gompertz") 
bouts <- readRDS("all.trajectory.daily.trips.DF.v10.dayevents.rds")


install.packages("dplyr")
library(dplyr)

table(bouts$id)
bouts$homerange.last<-0

animal<-"c001" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
plot(animal_day_subs$minsday,animal_day_subs$R2n,type="l")

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}



animal<-"c002" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
plot(animal_day_subs$minsday,animal_day_subs$R2n,type="l")

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c003" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c
for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c003.w" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c004" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c004.w" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c005" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c006" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c007" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c008" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c009" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c010" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c010.f" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c011" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c011.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c012" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c012.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c013" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c013.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c014" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c014.f" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}


animal<-"c015" #select an animal in the dataframe to include in the loop for plotting daily displacements 
# animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
plot(animal_day_subs$minsday,animal_day_subs$R2n)

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- last(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.last"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}


saveRDS(bouts, file = "all.trajectory.daily.trips.DF.v11.dayevents.rds")

#################################################################################
#####update number of observations
#################################################################################

###a loop to drop all the trips with observations less than required number 
bouts$nobs<-0

animal<-"c001" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c002" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}
animal<-"c003" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c003.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c004" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c004.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}
animal<-"c005" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c006" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c007" #select an animal in the dataframe to include in the loop for plotting daily displacements 
# animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
# plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c008" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c009" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c010" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c010.f" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c011" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c011.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c012" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c012.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c013" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c013.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c014" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c014.f" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c015" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

saveRDS(bouts, file = "all.trajectory.daily.trips.DF.v12.dayevents.rds")


hist(bouts$nobs) 
hist(bouts$nobs[bouts$nobs>0]) 
hist(bouts$nobs[bouts$nobs>0 & bouts$nobs<80])
#drop all the trips with < 80 observations 


#####################################################################
#### getting the first homerange 
#####################################################################
bouts <- readRDS("all.trajectory.daily.trips.DF.v12.dayevents.rds")
bouts$homerange.first<-0

animal<-"c001" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
plot(animal_day_subs$minsday,animal_day_subs$R2n,type="l")

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}



animal<-"c002" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
plot(animal_day_subs$minsday,animal_day_subs$R2n,type="l")

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c003" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c
for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c003.w" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c004" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c004.w" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c005" 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c006" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c007" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c008" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c009" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c010" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c010.f" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c011" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c011.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c012" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c012.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c013" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c013.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c014" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c014.f" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n,type="n", xlim = c(200,1100), ylim=c(0,200000000))c

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}


animal<-"c015" #select an animal in the dataframe to include in the loop for plotting daily displacements 
# animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
plot(animal_day_subs$minsday,animal_day_subs$R2n)

for (did in unique(bouts$dayid[bouts$id==animal])){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  print(sprintf("%s: %i %i observations",did,animal_day_subs$nobs[1],nrow(animal_day_subs)))
  if(nrow(animal_day_subs)>0)
  {
    m <- first(animal_day_subs$homerange)
  }
  #browser() 
  bouts[bouts$id==animal & bouts$dayid==did,"homerange.first"]<-m #nrow(animal_day_subs)
  #options(error=browser) 
  
  # if (coef(summary(m))[2,1]>0){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}


saveRDS(bouts, file = "all.trajectory.daily.trips.DF.v13.dayevents.rds")

###########################################################################
###recode to identify and classify data into permanent and temporary camps
#create a new field and assign it a unique id for camp types; 1. base camp and 2. temporary camp
#is the return at the permanent home or temprary camp? 
########################################################################

library(foreign)
library(lattice)
library(car)
require(rgdal)
require(raster) 
library(dplyr)

setwd("H:/Fern153/Gompertz") 

#homerange.last -the last points in the return trip and are within a 1000 m of home areas. 
#homerange.first-the start points in the outward trip and are within a 1000 me of the home range. 

bouts <- readRDS("all.trajectory.daily.trips.DF.v13.dayevents.rds")
bouts$camps <- "temporarily"
bouts$camps[bouts$homerange.last==1000 & bouts$homerange.first==1000] <- "permanent" # 


          #return1 return10 return11 return12 return13 return2 return3 return4 return5 return6 return7 return8 return9
#permanent 434975     2748      980      147        0  431869  471796  120335  128698  119417   25156   12098    3053
#temporarily123971    6740     8660     9243     9600  120106   91083   22665   16883   24203   11106   25837   34679

saveRDS(bouts, file = "all.trajectory.daily.trips.DF.v14.dayevents.rds")

##recoding variable label camps

bouts <- readRDS("all.trajectory.daily.trips.DF.v14.dayevents.rds")
table(bouts$camps)

bouts$camps[bouts$trips=="return12" & bouts$camps=="permanent"] <- "permastrange1" #no return trip & end in other homes
bouts$camps[bouts$trips=="return11" & bouts$camps=="permanent"] <- "permastrange1"  #no return trip & in other homes
bouts$camps[bouts$trips=="return10" & bouts$camps=="permanent"] <- "permastrange2"# a return trip but to neighbors 

bouts$camps[bouts$trips=="return13" & bouts$camps=="temporarily"] <- "tempostrange1" 
bouts$camps[bouts$trips=="return12" & bouts$camps=="temporarily"] <- "tempostrange2"
bouts$camps[bouts$trips=="return11" & bouts$camps=="temporarily"] <- "tempostrange3"
bouts$camps[bouts$trips=="return10" & bouts$camps=="temporarily"] <- "tempostrange4"


#bouts$camps[bouts$trips=="return2" & bouts$homerange==1000] <- 1 # 
##strange meaning gps left home either early or late but ended up in a permanent home elsewhere e.g. in the neighborhoods.a case on a day gps are collected from herders. 

#permanent camps 
test1 <- bouts[bouts$trips=="return12" & bouts$camps=="permanent",] 
test2 <- bouts[bouts$trips=="return11" & bouts$camps=="permanent",] 
test3 <- bouts[bouts$trips=="return10" & bouts$camps=="permanent",]#strange, trips returned to the closest neighbor. 
test4 <- bouts[bouts$trips=="return9" & bouts$camps=="permanent",] #fine though several trips have not gone far! 
test5 <- bouts[bouts$trips=="return8" & bouts$camps=="permanent",]#okay to keep 

####tempo camps 

test6 <- bouts[bouts$trips=="return13" & bouts$camps=="temporarily",]# 
test7 <- bouts[bouts$trips=="return12" & bouts$camps=="temporarily",]#
test8 <- bouts[bouts$trips=="return11" & bouts$camps=="temporarily",]#
test9 <- bouts[bouts$trips=="return10" & bouts$camps=="temporarily",]#
test10 <- bouts[bouts$trips=="return9" & bouts$camps=="temporarily",]#okay to keep

##permanent 
plot(test1$minsday,test1$R2n)
plot(test2$minsday,test2$R2n)
plot(test3$minsday,test3$R2n)
plot(test4$minsday,test4$R2n)
plot(test5$minsday,test5$R2n)

##tempo 
plot(test6$minsday,test6$R2n) #temp, not return to the start place 
plot(test7$minsday,test7$R2n)
plot(test8$minsday,test8$R2n)
plot(test9$minsday,test9$R2n)

plot(test10$minsday,test10$R2n) #a perfect returns to the start place  


table(bouts$camps)
xyplot(R2n~minsday|camps, data=bouts)
table(bouts$camps, bouts$trips)

saveRDS(bouts, file = "all.trajectory.daily.trips.DF.v15.dayevents.rds") 

##trancate data using nobs (number of rows in a burst) 
range(bouts$nobs)
hist(bouts$nobs) 
hist(bouts$nobs[bouts$nobs>0]) 
hist(bouts$nobs[bouts$nobs>0 & bouts$nobs<80])
#drop all the trips with < 80 observations 

#thresholding 
thresh <- 80

bouts.1 <- bouts[bouts$nobs>thresh,]
table(bouts.1$camps)
table(bouts.1$camps, bouts.1$trips)

xyplot(R2n~minsday|camps, data=bouts.1)

min(bouts.1$nobs)
hist(bouts.1$nobs) 
sum(is.na(bouts.1$R2n))
sum(is.na(bouts.1$nobs))

ta<- table(bouts.1$dayid,bouts.1$id)
write.csv(ta, file = "animal_codes_updated20180425_grt80nobs.csv")

head(bouts.1[[1]])
head(bouts.1[1])
# 
# #checks 
# length(bouts$x)
# length(bouts.1$x)
# length(bouts$x)-length(bouts.1$x)
# 
# bouts1<-bouts[bouts$NSD.max<thresh,]
#long.dry <- bouts.subset[bouts.subset$season=="Long Dry",]#subsetting

saveRDS(bouts.1, file = "all.trajectory.daily.trips.DF.v16.dayevents.rds") # 8 categories are presented and some may be dropped e.g. tempostrange1, tempostrange2, tempostrange3, and tempostrange4 and permastrange1, permastrange2. data shared with manuel on 25-04-2018 

###daily trips containing <80 observations


thresh <- 80
bouts.2 <- bouts[bouts$nobs<thresh,]

table(bouts.2$camps)
table(bouts.2$camps, bouts.2$trips)
xyplot(R2n~minsday|camps, data=bouts.2)

saveRDS(bouts.2, file = "all.trajectory.daily.trips.DF.v17.dayevents.rds") #obervations less than 80

bouts <- readRDS("all.trajectory.daily.trips.DF.v16.dayevents.rds")#dataset I shared with Manuel on 28.04.2018

##############################03-07-2018
###re-evaluating trips using different thresholds of nobs 
####03-07-2018
setwd("H:/Fern153/Gompertz/PastoGPS/29062018") 
bouts <- readRDS("bouts.day.180_1200.rds") #18 hours length for day events is a lot.


#################################################################################
#####update number of observations
#################################################################################
##data set: bouts.day. dataset exclude components on camp types
###a loop to drop all the trips with observations less than required number

#bouts.day <- readRDS("all.trajectory.daily.trips.DF.v9.b.dayevents.rds")
#bouts <- bouts.day

bouts$nobs<-0

animal<-"c001" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c002" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}
animal<-"c003" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c003.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c004" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c004.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}
animal<-"c005" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c006" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c007" #select an animal in the dataframe to include in the loop for plotting daily displacements 
# animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
# plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c008" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c009" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c010" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c010.f" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c011" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c011.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c012" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c012.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c013" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c013.w" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c014" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c014.f" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

animal<-"c015" #select an animal in the dataframe to include in the loop for plotting daily displacements 
animal_day_subs<-bouts[bouts$id==animal & bouts$dayid=="001_2012",]
#plot(animal_day_subs$minsday,animal_day_subs$R2n, xlim = c(0,1450), ylim=c(0,2629025910))

for (did in unique(bouts$dayid)){
  animal_day_subs<-bouts[bouts$id==animal & bouts$dayid==did,]
  bouts[bouts$id==animal & bouts$dayid==did,"nobs"]<-nrow(animal_day_subs)
  print(sprintf("%s: %i observations",did,nrow(animal_day_subs)))
  # if (nrow(animal_day_subs)>0 & nrow(animal_day_subs)<80){
  #   lines(animal_day_subs$minsday,animal_day_subs$R2n)
  # }
}

saveRDS(bouts, file = "all.trajectory.daily.trips.DF.v12.c.dayevents.rds")

##trancate data using nobs (number of rows in a burst) 
range(bouts$nobs)
hist(bouts$nobs) 
hist(bouts$nobs[bouts$nobs>0]) 
hist(bouts$nobs[bouts$nobs>0 & bouts$nobs<80])

fewerobs <- bouts$nobs[which(bouts$nobs>0 & bouts$nobs<80)]
summary(bouts$id[fewerobs])

thresh <- 140
p <- bouts[bouts$nobs<thresh,]
plot(p$minsday,p$R2n)

#drop all the trips with < 80 observations 

#thresholding 
thresh <- 60

bouts.1 <- bouts[bouts$nobs>thresh,]

min(bouts.1$nobs)
hist(bouts.1$nobs) 
hist(bouts$nobs[bouts$nobs>80 & bouts$nobs<140])
plot(bouts.1$minsday,bouts.1$R2n)

sum(is.na(bouts.1$R2n))
sum(is.na(bouts.1$nobs))

saveRDS(bouts.1, file = "all.trajectory.daily.trips.DF.v16.c.dayevents180to1200.rds") 
# data range include 1200<minsday>180. 

#bouts <- readRDS("all.trajectory.daily.trips.DF.v16.dayevents.rds")#dataset I shared with Manuel on 28.04.2018 

##
