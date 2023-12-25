#-----------------------------
# Plots
#-----------------------------

# INITS
library(sf)
library(rgdal) 
#setwd("~/mnt/agroscope_os/2/2/5/2/3/1935/Shibia")
setwd("H:/Fern153/Gompertz/PastoGPS/30082018") 
#Load functions
source("PastoGPS functions.R")
#try trip 5008
#Read-in data
d <- droplevels(readRDS("all.trajectory.daily.trips.DF.v16.b.dayevents.rds"))
m <- readGDAL("IntegrationV6_R4.tif")
res_all_df <- readRDS("res_all.rds")
obj <- res_all_df[6027,]# tested 6027, 5008
pt <- getSpatRef(obj)

  sub <- d[d$burst==as.character(obj$burst),]
  sub <- sub[!is.na(sub$x),]
  sub$minsday <- as.POSIXlt(sub$date)$hour*60+as.POSIXlt(sub$date)$min #Calculate mins per day
  pts <- st_as_sf(sub, coords = c("x", "y"))  
  sub$dis <- as.numeric(st_distance(pts[1,], pts)) #Calculate distance from 1st point
  if (!is.na(obj$asym)) sub$fDis <- with(obj, asym * (exp(-exp((b1-sub$minsday)/k1)) - exp(-exp((b2-sub$minsday)/k2))))

pdf("Figure2.pdf", 6,4)  
  layout(matrix(c(1,2), 1,2,byrow=F), width=c(1,1))
  par(mar=c(3,3,1,1), mgp=c(2,.5,0))
  lns <- st_linestring(as.matrix(sub[c("x", "y")]))
  m_y<-mean(st_bbox(lns)[c(2,4)])
  plot(lns, ylim=m_y+c(3000,-3000), axes=T, xlab="X (m)", ylab="Y (m)")
# image(m, add=T, col=rainbow(10))
# plot(lns, add=T)
  plot(st_geometry(pts[c(1,nrow(pts)),]), add=T, col="black", pch=c(16,16))
  plot(st_geometry(pt[3:4,]), col="red", pch=16, add=T)
  lines(st_coordinates(pt[3:4,]),col="red")
  title("a", adj=0.05, line = -2) 

  plot(sub$date, sub$dis , t="l", ylab="Distance from start (m)", xlab="Time")
  lines(sub$date, sub$fDis, col="blue", lwd=3)
  t1 <- as.POSIXct(obj$GrazStart*60, origin=strftime(sub$date[1], format="%Y-%m-%d"))
  x1 <-  with(obj, asym * (exp(-exp((b1-obj$GrazStart)/k1)) - exp(-exp((b2-obj$GrazStart)/k2))))
  arrows(t1,x1,t1,0, length=.1, col="red")
  t2 <- as.POSIXct(obj$GrazEnd*60, origin=strftime(sub$date[1], format="%Y-%m-%d"))
  x2 <-  with(obj, asym * (exp(-exp((b1-obj$GrazEnd)/k1)) - exp(-exp((b2-obj$GrazEnd)/k2))))
  arrows(t2,x2,t2,0, length=.1, col="red")
  title("b", adj=0.05, line = -2) 
dev.off()
