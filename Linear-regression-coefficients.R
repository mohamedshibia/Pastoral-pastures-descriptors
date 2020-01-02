#code supplimentary
#extract coefficients of linear regression analysis

setwd("") # set working directory
rm(list = ls()) 
 
require("rgdal") #install required packages 

ds <- GDAL.open("input_PVI_files_stacked") #import stacked PVIs from working directory
projection <- getProjectionRef(ds) 
gt <- .Call("RGDAL_GetGeoTransform", ds, PACKAGE = "rgdal") 

image <- getRasterData(ds) 
class(image) 
dim(image) 

nr <- nrow(image) 
nc <- ncol(image) 

slope <- matrix(nrow=nr, ncol=nc) #create empty matrix to the size of the images
p.val <- matrix(nrow=nr, ncol=nc)

nd <- 0
x <- c(0,16,69) # day counts start with the image on 20-12-2014, 05-01-2015, and 10-03-2015
for (r in 1:nr) {
  print(r)
  for (c in 1:nc) {
    y <- image[r,c,]
    if (is.na(y[1])){
      slope[r,c] <- nd
      p.val[r,c] <- nd
          }else{
    model<-lm(y~x) 
    s<-summary(model)
    slope[r,c]<-s$coefficients[2,1] 
    p.val[r,c]<-s$coefficients[2,4] 
    }
  }
}


driver <- new("GDALDriver", "GTiff") # pick preferred driver, and here we select geotiff one. 
ds.out <- new("GDALTransientDataset", driver, nc, nr, 2, type="Float32") 
putRasterData(ds.out, slope, band=1) 
putRasterData(ds.out, p.val, band=2) 
.Call("RGDAL_SetProject", ds.out, projection, PACKAGE = "rgdal") 
.Call("RGDAL_SetGeoTransform", ds.out, gt[1:6], PACKAGE = "rgdal")
saveDataset(ds.out, "Output.tif") #write slope and p.value regression coefficients.

GDAL.close(ds.out) 
GDAL.close(ds)
