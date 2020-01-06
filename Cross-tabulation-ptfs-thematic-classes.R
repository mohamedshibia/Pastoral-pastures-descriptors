## contingency table and chi- square statistics   

setwd()# set working directory 

#cross tabulation of map classes and plant functional types
tab <-matrix(c(71,31,11,67,14,3,0,0,0,0,14,8,0,7,0,0,0,1,42,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,35,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,4,0,15,15,16,46,7,2,1,1,0,5,20,0,2,0,0,0,2,1,3,0,20,14,5,14,64,0,1,27,4,35,11,0,1,0,8,6,78,5,7,16,0,0,30,4,9,4,0,0,0,0,0,0,0,0,37,0,0,2,8,0,5,0,0,0,0,0,0,0,0,29,0,38,24,0,7,19,9,7,0,1,153,25,0,7,0,0,0,0,0,10,0,0,0,0,0,0,0,0,23,166,0,0,0,0,0,19,11,0,23,0,0), ncol=11)

rownames(tab) <- c(111,121,122,131,132,133,134,143,144,145,231,232,245,332,345,445) # map class labels
colnames(tab) <- c("ba","dds","fdd","fde","gs","gt","sdd","sde","sdhd","sds","sts") # plant functional types (pft)
print(tab)

##column names (pft); ba:annuals, dds:  dwarf deciduous shrubs, fdd:  dense deciduous forest: fde: dense evergreen forests, gs: dense short grasses,gt: dense tall grasses, sdd: dense deciduous shrubs, sde: dense evergreen shrubs, sdhd: highly deciduous sparse shrubs, sds: sparse dwarf shrubs, sts: sparse tall shrubs 

#to determine strength of the links between pairs of plant functional types (pft) and thematic classes, calculate chi-square statistics and then extract standardized residuals. Cells with highest values contribute the most to the total chi-square. Positive ones denote attraction and negative values represent repulsion.
# 
chisq <- chisq.test(tab,simulate.p.value = TRUE, B = 9999) 

# Visualize cell chi-square residuals   

m<-as.matrix(chisq$residuals)
library(fields)
lev <- seq(-35,35,2)
colors <- c(colorRampPalette(c("dark red", "white", "dark blue"))(length(lev)-1))
image.plot(t(m)[,16:1], breaks = lev, col=colors, axes=F)
axis(2, rownames(m)[16:1], at=seq(0,1, length.out=16), tick=F, las=2)
axis(3, colnames(m), at=seq(0,1, length.out=11), tick=F)
