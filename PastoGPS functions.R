# FUNCTIONS
#----------

fitNLFun <- function(sub, model="Gompertz", nsd=F, aug=20, plot=F){
  # Input: sub is data of an individual trip, needs to be a data.frame with columns x, y, #date
  #       model is either "Gompertz", "Exponential" or "Exponential Upper"
  #       nsd is TRUE/FALSE if nsd is used or not
  #       plot is TRUE/FALSE if plot is generated or not
  # Output is named vector of estimated values
  type <- "NA coordinates"
  mig <- FALSE

  # 1. Transform input data into sf object
  if(length(which(!is.na(sub$x)))>20){ #Check that coordinates are available
  sub <- sub[!is.na(sub$x),]
  pts <- st_as_sf(sub, coords = c("x", "y"))
  sub$dis <- as.numeric(st_distance(pts[1,], pts)) #Calculate distance from 1st point
  sub$minsday <- as.POSIXlt(sub$date)$hour*60+as.POSIXlt(sub$date)$min #Calculate mins per day

  # Routine for cases were start camp and end camp are not identical 
  # Calculates the distance from a linear displacement between start and end camp
  if(rev(sub$dis)[1] > 1000){
    mig <- TRUE
    l <- nrow(sub)
    interp <- data.frame(minsday=seq(sub$minsday[1], rev(sub$minsday)[l],5))
    interp$x <- seq(sub$x[1], sub$x[l],length.out=nrow(interp))
    interp$y <- seq(sub$y[1], sub$y[l],length.out=nrow(interp))
    sub$dis <- sqrt((sub$x -interp$x[match(sub$minsday, interp$minsday)])^2 + (sub$y -interp$y[match(sub$minsday, interp$minsday)])^2)
  }
  
  # 2. Fit the model (try is used for cases of non conversion)
  if(model=="Exponential"){ #Double exponential
    if(nsd) M <- try(nls(dis^2 ~  asym /(1+exp((b1-minsday)/k1)) + (-asym /(1 + exp((b2-minsday)/k2))), data=sub, start = c(asym=5000^2,b1=300,b2=800,k1=50,k2=70)), silent=T) #Try is used to continue in case of non-convergence
    else M <- try(nls(dis ~  asym /(1+exp((b1-minsday)/k1)) + (-asym /(1 + exp((b2-minsday)/k2))), data=sub, start = c(asym=5000,b1=300,b2=800,k1=50,k2=70)), silent=T) #Try is used to continue in case of non-convergence
  }
  if(model == "Exponential Upper"){ #Double exponential with upper bound
    if(nsd) M <- try(nls(dis^2 ~  asym /(1+exp((b1-minsday)/k1)) + (-asym /(1 + exp((b2-minsday)/k2))), data=sub, start = c(asym=5000^2,b1=300,b2=800,k1=50,k2=70), upper=c(asym=max(sub$dis)), algorithm="port"), silent=T)
    else M <- try(nls(dis ~  asym /(1+exp((b1-minsday)/k1)) + (-asym /(1 + exp((b2-minsday)/k2))), data=sub, start = c(asym=5000,b1=300,b2=800,k1=50,k2=70), upper=c(asym=max(sub$dis)), algorithm="port"), silent=T)
  }
  if(model == "Gompertz"){ #Double Gompertz
    if(nsd) M <- try(nls(dis^2 ~   asym * (exp(-exp((b1-minsday)/k1)) - exp(-exp((b2-minsday)/k2))), data=sub, start = c(asym=5000^2,b1=300,b2=800,k1=80,k2=5), control=list(maxiter=5000,minFactor=1/1000000)), silent=T)
    else{ M <- try(nls(dis ~   asym * (exp(-exp((b1-minsday)/k1)) - exp(-exp((b2-minsday)/k2))), data=sub, start = c(asym=5000,b1=300,b2=800,k1=80,k2=5), control=list(maxiter=5000,minFactor=1/1000000)), silent=T)
    # Can have a series of calls with different starting values here for those cases,in which the default does not fit
    if(class(M) == "try-error") M <- try(nls(dis ~   asym * (exp(-exp((b1-minsday)/k1)) - exp(-exp((b2-minsday)/k2))), data=sub, start = c(asym=1000,b1=10,b2=800,k1=80,k2=5), control=list(maxiter=5000,minFactor=1/1000000)), silent=T)
    if(class(M) == "try-error") M <- try(nls(dis ~   asym * (exp(-exp((b1-minsday)/k1)) - exp(-exp((b2-minsday)/k2))), data=sub, start = c(asym=5000,b1=300,b2=600,k1=5,k2=5), control=list(maxiter=5000,minFactor=1/1000000)), silent=T)
    if(class(M) == "try-error") M <- try(nls(dis ~   asym * (exp(-exp((b1-minsday)/k1)) - exp(-exp((b2-minsday)/k2))), data=sub, start = c(asym=10000,b1=10,b2=1000,k1=50,k2=30), control=list(maxiter=5000,minFactor=1/1000000)), silent=T)
    if(class(M) == "try-error"){
      # Augment data before and after start
      sub2 <- sub[,c("dis", "minsday")]
      sub2 <- rbind(data.frame(dis=rep(0, aug), minsday=seq(sub$minsday[1], by=-5, length=aug)),
                    sub2,
                    data.frame(dis=rep(sub$dis[nrow(sub)], aug), minsday=seq(sub$minsday[nrow(sub)], by=5, length=aug)))
      M <- try(nls(dis ~   asym * (exp(-exp((b1-minsday)/k1)) - exp(-exp((b2-minsday)/k2))), data=sub2, start = c(asym=5000,b1=300,b2=800,k1=80,k2=5), control=list(maxiter=5000,minFactor=1/1000000)), silent=T)
    }  
    }
  }
  
  # 3. Extract parameters
  if(class(M) == "try-error") M_co <- rep(NA,5) else M_co <- coef(M)
  names(M_co) <- c("asym", "b1", "b2", "k1", "k2")
  
  # standard error of the estimate for population data
  # see http://r.789695.n4.nabble.com/How-to-calculate-standard-error-of-estimate-S-for-my-non-linear-regression-model-td4712808.html
  if(class(M) == "try-error") sigma.est <- NA else sigma.est <- sqrt(sum(resid(M)^2) / nrow(sub))
  if(class(M) == "try-error") type<-"try-error" else type <- "OK"
  
  out <- cbind(data.frame(id = sub$id[1], burst=sub$burst[1], date=format(sub$date[1], format="%d.%m.%y"),
                          start_t=format(sub$date[1], format="%H:%M"), end_t=format(sub$date[nrow(sub)], format="%H:%M"),
                          start_x = sub$x[1], start_y=sub$y[1], end_x= sub$x[length(sub$x)], end_y = sub$y[length(sub$y)],
                          model=model, nsd=nsd, sigma.est=sigma.est, migratory=mig, type=type),
               data.frame(t(M_co)))
  } else out <- data.frame(id = NA, burst=NA, date=NA, start_t=NA, end_t=NA,
                         start_x = NA, start_y=NA, end_x= NA, end_y=NA,
                         model=NA, nsd=NA, sigma.est=NA, migratory=NA, type=type, asym=NA, b1=NA, b2=NA, k1=NA, k2=NA)

  # 4. Possibly do a plot
  if(plot){
    layout(matrix(c(1,1,2,3), 2,2,byrow=F), width=c(2,1))
    lns <- st_linestring(as.matrix(sub[c("x", "y")]))
    m_y<-mean(st_bbox(lns)[c(2,4)])
    plot(lns, ylim=m_y+c(3000,-3000), axes=T)
    title(paste(sub$burst[1],":", sub$Trip[1], ":", format(sub$date[1], format="%d.%m.%y")))
    if(class(M) != "try-error") legend("topleft", legend=round(coef(M),0))
    if(class(M) != "try-error") ymax <- abs(M_co[1]) else ymax <- 3000
    plot(sub$date, sub$dis , t="l",ylim=c(0, ymax), ylab="Distance from start")
    if (class(M) != "try-error") sub$fDis <- with(data.frame(t(coef(M))), asym * (exp(-exp((b1-sub$minsday)/k1)) - exp(-exp((b2-sub$minsday)/k2))))
    if(class(M) != "try-error") lines(sub$date, sub$fDis, col="red", lwd=2)
    plot(sub$date, c(0,cumsum(st_distance(pts$geometry[-nrow(pts)], pts$geometry[-1], by_element=T))), t="l", ylim=c(0,20000), ylab="Commulated distance")
  }
  
  return(out)
}


evaluateTrip <- function(obj){
  # 1. Predict values
  minsday <- seq(0, 24*60, 1)
  if (!is.na(obj$asym)) {fDis <- with(obj, asym * (exp(-exp((b1-minsday)/k1)) - exp(-exp((b2-minsday)/k2))))
  
  # 2. Key points
  # Trip starts and ends if distance is >100m from initial point
  TripStart <- minsday[which(fDis>100)[1]]
  TripEnd <- minsday[rev(which(fDis>100))[1]]
  
  # Grazing starts and ends if distance is >100m from maximum
  maxDis <- max(fDis)
  GrazStart <- minsday[which(maxDis-fDis<100)[1]]
  GrazEnd <- minsday[rev(which(maxDis-fDis<100))[1]]
  # 3. Derived quantities 
  TripLength <- TripEnd - TripStart
  GrazLength <- GrazEnd - GrazStart
  return(cbind(obj, data.frame(TripStart, TripEnd, GrazStart, GrazEnd, maxDis, TripLength, GrazLength)))
  } else return(cbind(obj, data.frame(TripStart=NA, TripEnd=NA, GrazStart=NA, GrazEnd=NA, maxDis=NA, TripLength=NA, GrazLength=NA)))
}

# Getting spatial reference for key points
getSpatRef <- function(obj, plot=F){
  if(obj$type != "NA coordinates"){
  sub <- d[d$burst==as.character(obj$burst),]
  sub <- sub[!is.na(sub$x),]
  sub$minsday <- as.POSIXlt(sub$date)$hour*60+as.POSIXlt(sub$date)$min #Calculate mins per day
  pts <- st_as_sf(sub, coords = c("x", "y"))  
  if(!is.na(obj$TripStart)){ keyPts <- pts[c(which(pts$minsday>obj$TripStart)[1], rev(which(pts$minsday<obj$TripEnd))[1],
                                             which(pts$minsday>obj$GrazStart)[1], rev(which(pts$minsday<obj$GrazEnd))[1]),]
  keyPts$Name <- c("TripStart", "TripEnd", "GrazStart", "GrazEnd")
  } else keyPts=NA
  if(anyNA(keyPts)) keyPts=NA
  if(!is.na(obj$TripStart) && plot==T){
    lns <- st_linestring(as.matrix(sub[c("x", "y")]))
    m_y<-mean(st_bbox(lns)[c(2,4)])
    plot(lns, ylim=m_y+c(3000,-3000), axes=T)
    plot(st_geometry(keyPts), add=T, pch=c(16,16,17,17), col=1:2, cex=2)
  }
  } else keyPts=NA
  return(keyPts)
}

# Convert key points to line
points2line <- function(obj){
  if(class(obj)[1]=="sf") line <- st_sf(st_cast(st_combine(obj[3:4,]), "LINESTRING")) else line <- NA
  if(!is.na(line)){
    line$id <- obj[3,]$id
    line$burst <- obj[3,]$burst
    line$dateStart <- obj[3,]$date
    line$dateEnd <- obj[3,]$date
  }
  return(line)
}





plotTrip <- function(obj, pts=NULL){
  sub <- d[d$burst==as.character(obj$burst),]
  sub <- sub[!is.na(sub$x),]
  sub$minsday <- as.POSIXlt(sub$date)$hour*60+as.POSIXlt(sub$date)$min #Calculate mins per day
  pts <- st_as_sf(sub, coords = c("x", "y"))  
  sub$dis <- as.numeric(st_distance(pts[1,], pts)) #Calculate distance from 1st point
  if (!is.na(obj$asym)) sub$fDis <- with(obj, asym * (exp(-exp((b1-sub$minsday)/k1)) - exp(-exp((b2-sub$minsday)/k2))))
  
  if(rev(sub$dis)[1] > 1000){
    l <- nrow(sub)
    interp <- data.frame(minsday=seq(sub$minsday[1], sub$minsday[l],5))
    interp$x <- seq(sub$x[1], sub$x[l],length.out=nrow(interp))
    interp$y <- seq(sub$y[1], sub$y[l],length.out=nrow(interp))
    sub$dis2 <- sqrt((sub$x -interp$x[match(sub$minsday, interp$minsday)])^2 + (sub$y -interp$y[match(sub$minsday, interp$minsday)])^2)
  }
  
  layout(matrix(c(1,1,2,3), 2,2,byrow=F), width=c(2,1))
  lns <- st_linestring(as.matrix(sub[c("x", "y")]))
  m_y<-mean(st_bbox(lns)[c(2,4)])
  plot(lns, ylim=m_y+c(3000,-3000), axes=T)
  plot(st_geometry(pts[c(1,nrow(pts)),]), add=T, col="red", pch=c(16,17))
  title(paste(sub$burst[1],":", sub$Trip[1], ":", format(sub$date[1], format="%d.%m.%y")))
  if(!is.na(obj$asym)) legend("topleft", legend=round(obj[12:16],0))
  if(!is.na(obj$asym)) ymax <- abs(obj$asym) else ymax <- max(sub$dis)
  par(mar=c(2,2,1,1))
  plot(sub$date, sub$dis , t="l",ylim=c(0, ymax), ylab="Distance from start")
  if(rev(sub$dis)[1] > 1000) lines(sub$date, sub$dis2, t="l",lwd=2, lty=2)
  if(!is.na(obj$asym)) lines(sub$date, sub$fDis, col="red", lwd=2)
  plot(sub$date, c(0,cumsum(st_distance(pts$geometry[-nrow(pts)], pts$geometry[-1], by_element=T))), t="l", ylim=c(0,20000), ylab="Commulated distance")
}
