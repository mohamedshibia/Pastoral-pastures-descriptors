library(compositions)
source("PastoGPS functions.R")

res_all_df <- readRDS("res_all_df_updated.rds")
res_all_df$year <- as.POSIXlt(res_all_df$date)$year
res_all_df$month <- as.POSIXlt(res_all_df$date)$mon

#Calculate compostions per month
group <- factor(paste(res_all_df$id, res_all_df$year, res_all_df$month, sep="_"))#674 combinations 
tab <- do.call("rbind",tapply(factor(res_all_df$pam.PC.cluster), group, function(x) as.numeric(table(x))))
table(rowSums(tab)) # I notice that there are many months with only a few days recorded.

Y <- acomp(tab/rowSums(tab))
gr <- data.frame(group = levels(group))
gr$id <- sapply(gr$group, function(x) strsplit(as.character(x), "_")[[1]][1])
gr$year <- sapply(gr$group, function(x) strsplit(as.character(x), "_")[[1]][2])
gr$month <- sapply(gr$group, function(x) strsplit(as.character(x), "_")[[1]][3])
gr$site <- res_all_df$site[match(gr$id, res_all_df$id)]
gr$herd <- res_all_df$herd[match(gr$id, res_all_df$id)]
gr$TLU <- as.numeric(as.character(res_all_df$TLU[match(gr$id, res_all_df$id)]))
gr$season <- res_all_df$season[match(gr$month, res_all_df$month)]

M <- lm(ilr(Y) ~ gr$TLU + gr$site + gr$season)
summary(M)
anova(M)
ilrInv(coef(M), orig=Y)
