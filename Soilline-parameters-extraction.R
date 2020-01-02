# Extract Soilline and soilline parameters
# R-Skript: 
rm(list=ls())
library(rgdal) # install required packages


setwd("") # set working directory 
dir()
list<-dir()


# adapt accordingly
list_num <- 21 #select the file position number 21 in the directory e.g. "LC81680572014354.dat" 

quant <- 0.05  #percent quantile of pixel to be included in the analysis. One can increase. Increasing quantile-th will subsequently increase computation time to complete processing. 
airange <- seq(500,2500,by=100) #range in the sequence; in our case it'll start at 500 and end with 2500 of RED values and a sequence of 100.


#read data, clean data (remove 0 values) and sort data
x <- readGDAL(list[list_num])$band3 #Red
y <- readGDAL(list[list_num])$band4 #NIR

ysub <- subset(y, y > 0 & x > 0 )
xsub <- subset(x, y > 0 & x > 0 )

sortx <- sort(xsub, index.return=TRUE)
sorty <- ysub[sortx[[2]]]

min(xsub)
max(xsub)

# adapt boundaries and intervals according to value range

length(airange)

xsortvalue <- sortx[[1]]

xmedian <- c()
y005 <- c()

# calculation of quantiles

for (i in 1:(length(airange)-1)) {
	xmedian[i] <-median(xsortvalue[(xsortvalue > airange[i] & xsortvalue < airange[i+1])], na.rm = TRUE)
	y005[i] <- quantile(sorty[(xsortvalue > airange[i] & xsortvalue < airange[i+1])], quant,na.rm = TRUE)

}

# fit soil line
lowerfit = nls(y005 ~ a + (b*xmedian), start = list(a=0, b=1))
coefficients(lowerfit)
summary(lowerfit)

# plot: name of plot corresponds to image name
output_plot <- paste("Plot_", unlist(strsplit(list[list_num], ".",fixed=TRUE))[1], ".png", sep="")
png(filename= output_plot, width = 400, height = 400)

plot(xsub, ysub, xlim=c(0,4000), ylim=c(0,5000), main="(III)", xlab = "Red", ylab="NIR")
points(xmedian, y005, col="red")

xmedianNA <- xmedian[(!is.na(xmedian)) & xmedian>0]
lines(xmedianNA,predict(lowerfit), col = "red")##get some estimation of goodness of fit


mtext(paste("a: ", round(as.numeric(coefficients(lowerfit)[1]), 3), sep=""),3,2, adj=0.05)
mtext(paste("b: ", round(as.numeric(coefficients(lowerfit)[2]), 3), sep=""),3,1, adj=0.05)


dev.off()

# # print coefficients
coefficients(lowerfit)
