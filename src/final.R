# Working directory
switch(Sys.info()[["sysname"]], 
       "Linux" = setwd("D:/active/kilimanjaro"), 
       "Windows" = setwd("D:/active/kilimanjaro"))

library(kza)
library(latticeExtra)
library(Kendall)
library(foreach)
library(ade4)
library(reshape2)
library(ggplot2)
source("vectorHarmonics.R")

outLayer <- function(x, y) {
  x + as.layer(y)
}


#### Long-term precipitation analysis ##########################################
# Compute long-term anomalies and trends for KIA precipitation and create
# publication quality figure.

# Read precipitation data and create a continous time series.
precip <- 
  read.csv("kilimanjaro_gsod_dynamics/gsod_precip/kia_prcp_1975_2014_mnthly.csv", 
           stringsAsFactors = FALSE)
precip[, 1] <- as.Date(precip[, 1])
st <- precip[1, 1]
nd <- precip[nrow(precip), 1]
ts <- seq(st, nd, "month")
precip <- merge(data.frame(ts), precip, by = 1, all.x = TRUE)

# Compute 3 month running mean of original precipitation values using a
# Kolmogorov-Zurbenko filter with one iteration
precip$kz03k01 <- kz(precip$P_RT_NRT, m = 3, k = 1)

# Compute deseasoned precipitation time series and corresponding 3 month running
# mean using a Kolmogorov-Zurbenko filter with two iterations to close gaps
precip$ssn <- precip$P_RT_NRT - rep(sapply(1:12, function(i) {
  mean(precip$P_RT_NRT[seq(i, nrow(precip), 12)], na.rm = TRUE)
}), nrow(precip) / 12)

precip$ssn_kz03k01 <- kz(precip$ssn, m = 3, k = 1)
precip$ssn_kz03k02 <- kz(precip$ssn, m = 3, k = 2)

clr <- as.character(ifelse(precip$ssn_kz03k02 > 0, "blue", "red"))
clr.blue <- brewer.pal(9, "Blues")[2]
clr.red <- brewer.pal(9, "Reds")[1]

# Create publication quality figure of long term precipitation trends
plot.precip.ssn_kz03k02 <- 
  xyplot(ssn_kz03k02 ~ ts, data = precip, origin = 0, type = "h",
         border = "transparent", col = clr, asp = 0.25,
         xlab = "", ylab = "Precipitation [mm]", 
         lwd = 1.5, ylim = c(-250, 350), as.table = TRUE,
         scales = list(x = list(axs = "i")),
         xscale.components = xscale.components.subticks,
         yscale.components = yscale.components.subticks,
         panel = function(x, y, ...) {
           panel.xblocks(x, y > 0, 
                         col = clr.blue)
           panel.xblocks(x, y < 0, 
                         col = clr.red)
           panel.xyplot(x, y, ...)
           panel.smoother(x, y, method = "lm", 
                          col = "black", 
                          col.se = "black",
                          alpha.se = 0.3, lty = 2)
         })

MannKendall(precip$ssn_kz03k02)



#### Seasonal precipitation analysis ###########################################
# Compute seasonal precipitation dynamics for each year of the time series
# centered arround Dezember and create publication quality figure. 

# Compute a three year running mean of the seasonal precipitation dynamics
precip.seasonal03yr <- 
  foreach(st = 1975:2003, nd = 1985:2013, .combine = "rbind") %do% {
    st.process <- grep(st, precip[, 1])[7]
    nd.process <- grep(nd, precip[, 1])[length(grep(nd, precip[, 1]))-6]
    precip.process <- precip[st.process:nd.process, ]
    data.frame(season = paste(st, nd, sep = "-"), 
               p_dyn = vectorHarmonics(precip.process$kz03k01, 
                                       st = c(st, 1), 
                                       nd = c(nd, 12), 
                                       m = 5))
  }

# Split the running mean data set by year and create publication quality figure;
# to smooth the cycle while not reducing the rainfall amounts significantly, 
# use a spline prediction.
precip.seasonal03yr.split <- split(precip.seasonal03yr, precip.seasonal03yr$season)
colors <- colorRampPalette(brewer.pal(11, "RdYlGn"))

xhres <- seq(1, 12, 0.01)
at <- seq(1, 1101, 100)
labels <- c(seq(7, 12, 1), seq(1, 6, 1))
plot.precip.seasonal03yr.split <- 
  lapply(seq(precip.seasonal03yr.split), function(x){
    sp <- predict(smooth.spline(
      precip.seasonal03yr.split[[x]]$p_dyn, spar=0.01), xhres)
    xnew <- as.factor(sp[[1]])
    factor(xnew, levels = c(seq(7,12,0.01), seq(1,6,0.01)))
    xyplot(sp[[2]] ~ xnew, type = "l", lwd = "2", 
           ylim = c(0, 170), 
           col = colors(length(precip.seasonal03yr.split))[x],
           scale=list(x=list(at = at, labels = labels)),
           xlab = "Month", ylab = "Precipitation")
  })
plot.precip.seasonal03yr.all <- 
  Reduce("outLayer", plot.precip.seasonal03yr.split)

# Non-spline version
# plot.precip.seasonal03yr.split <- 
#   lapply(seq(precip.seasonal03yr.split), function(i) {
#     xyplot(p_dyn ~ factor(c(7:12, 1:6), levels = c(7:12, 1:6)), 
#            data = precip.seasonal03yr.split[[i]], type=c("l"), 
#            col = colors(length(precip.seasonal03yr.split))[i], 
#            ylim = c(0,170),
#            xlab = "Month", ylab = "Precipitation")
#   })
# plot.precip.seasonal03yr.all <- 
#   Reduce("outLayer", plot.precip.seasonal03yr.split)



#### Precipitation analysis vs ENSO ############################################
# Prepare ONI record and classifiy years as La Nina (L), El Nino (E) or 
# normal (N); weak ENSO cycles are classified as normal
oni <- read.csv("enso/enso_and_iod.csv", skip = 1, header = TRUE)
oni$Season <- paste0(oni$Season, oni$X, oni$X.1)
oni <- oni[, -grep("X", names(oni))]

oni$TypeClass <- "N"
oni$TypeClass[grep("L", oni$Type)] <- "LS"
oni$TypeClass[grep("WL", oni$Type)] <- "LW"
oni$TypeClass[grep("E", oni$Type)] <- "ES"
oni$TypeClass[grep("WE", oni$Type)] <- "EW"
# oni$TypeClass[grep("W", oni$Type)] <- "N"
# oni$TypeClass[grep("M", oni$Type)] <- "X"

oni.reshape <- do.call("rbind", lapply(seq_len(nrow(oni)), function(i) {
  st <- as.Date(paste0(substr(oni[i, 3], 1, 4), "-07-01"))
  nd <- as.Date(paste0(substr(oni[i, 3], 6, 9), "-06-01"))
  data.frame(oni[i, c(1:3, ncol(oni))], 
             oni = as.numeric(oni[i, 4:(ncol(oni)-1)]), 
             month = seq(st, nd, "month"))
}))

# Adjust precipitation data records to ENSO cycle logic (i.e. start in July and
# end in June) and combine precipitation and ONI data set
precip.shift06m <- precip[7:(nrow(precip)-6), ]

precip.shift06m.oni <- 
  merge(precip.shift06m, oni.reshape, all.x = TRUE, by.x = "ts", by.y = "month")

precip.shift06m.oni.split <- 
  split(precip.shift06m.oni, precip.shift06m.oni$TypeClass)

precip.shift06m.oni.split.median <- 
  foreach(i = precip.shift06m.oni.split) %do% {
    sapply(1:12, function(j) {
      median(i$P_RT_NRT[seq(j, nrow(i), 12)], na.rm = TRUE)
    })
  }

# Produce publication quality graph for mean seasonal cycles during El Nino,
# La Nina and normal situations. To smooth the cycle while not reducing the
# rainfall amounts significantly, use a spline prediction.
if(length(precip.shift06m.oni.split) == 5){
  colors <- c("blue", "lightblue", "red", "bisque", "black")
} else {
  colors <- c("blue", "red", "black")  
}
xhres <- seq(1, 12, 0.01)
at <- seq(1, 1101, 100)
labels <- c(seq(7, 12, 1), seq(1, 6, 1))
plot.precip.shift06m.oni.split.median <- 
  lapply(seq(precip.shift06m.oni.split.median), function(x){
    sp <- predict(smooth.spline(
      precip.shift06m.oni.split.median[[x]], spar=0.01), xhres)
    xnew <- as.factor(sp[[1]])
    factor(xnew, levels = c(seq(7,12,0.01), seq(1,6,0.01)))
    xyplot(sp[[2]] ~ xnew, type = "l", lwd = "2", 
           ylim = c(0, 200), col = colors[x],
           scale=list(x=list(at = at, labels = labels)),
           xlab = "Month", ylab = "Precipitation")
  })
plot.precip.shift06m.oni.split.median.all <- 
  Reduce("outLayer", plot.precip.shift06m.oni.split.median)

# Non-spline version
# plot.precip.shift06m.oni.split.median <- 
#   lapply(seq(precip.shift06m.oni.split.median), function(x){
#     xyplot(precip.shift06m.oni.split.median[[x]] ~ 
#              factor(c(7:12, 1:6), levels = c(7:12, 1:6)), type = "l", 
#            ylim = c(0, 200), col = colors[x])
#   })
# 
# plot.precip.shift06m.oni.split.median.all <- 
#   Reduce("outLayer", plot.precip.shift06m.oni.split.median)
# plot.precip.shift06m.oni.split.median.all <- 
#   Reduce("outLayer", plot.precip.shift06m.oni.split.median)

# Just for background information: the individual seasons which form the above
# combined El Nino, La Nina and normal situation spline.
precip.shift06m.oni.split.season <- lapply(precip.shift06m.oni.split, function(x){
  split(x, x$Season)
})
xhres <- seq(1, 12, 0.01)
at <- seq(1, 1101, 100)
labels <- c(seq(7, 12, 1), seq(1, 6, 1))
plot.precip.shift06m.oni.split.season <- 
  lapply(seq(precip.shift06m.oni.split.season), function(x){
    lapply(seq(precip.shift06m.oni.split.season[[x]]), function(i) {
      sp <- with(precip.shift06m.oni.split.season[[x]][[i]][!is.na(precip.shift06m.oni.split.season[[x]][[i]]$P_RT_NRT),], smooth.spline(P_RT_NRT, spar=0.01))
      sp <- predict(sp, xhres)
    xnew <- as.factor(sp[[1]])
    factor(xnew, levels = c(seq(7,12,0.01), seq(1,6,0.01)))
    xyplot(sp[[2]] ~ xnew, type = "l", lwd = "2", 
           ylim = c(0, 700), col = colors[x],
           scale=list(x=list(at = at, labels = labels)),
           xlab = "Month", ylab = "Precipitation")
  })})
plot.precip.shift06m.oni.split.season.all <- 
  Reduce("outLayer", lapply(plot.precip.shift06m.oni.split.season, function(x){
    Reduce("outLayer", x)
  }))


# Alternative boxplot visualization
## bwplots
precip.shift06m.oni$month <- substr(precip.shift06m.oni$ts, 6, 7)
precip.shift06m.oni.melt <- melt(precip.shift06m.oni, 
                                 id.vars = c("month", "TypeClass"),
                                 measure.vars = "P_RT_NRT")
bw.plot <- ggplot(data = precip.shift06m.oni.melt, 
                  aes(y = value, x = month, fill = TypeClass, 
                      dodge = TypeClass))
bw.plot + geom_boxplot()

# Create publication quality figures of correlations between ONI and 
# precipitation
years <- unique(substr(precip$ts, 1, 4))
precip.shift06m.mat.ssn_kz03k01 <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.oni[grep(i, precip.shift06m.oni$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$ssn_kz03k01[j]
  }
}

oni.reshape.shift06m.mat <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.oni[grep(i, precip.shift06m.oni$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$oni[j]
  }
}

precip.shift06m.ssn_kz03k01.oni.reshape.shift06m.cor <- 
  cor(oni.reshape.shift06m.mat, precip.shift06m.mat.ssn_kz03k01, use = "complete.obs", 
      method = "kendall")
for (r in 1:12) {
  for (c in 1:12) {
    if (r > c) {
      precip.shift06m.ssn_kz03k01.oni.reshape.shift06m.cor[r,c] <- NA
    }
  }
}

colors <- colorRampPalette(brewer.pal(9, "RdBu"))
levelplot(t(precip.shift06m.ssn_kz03k01.oni.reshape.shift06m.cor), col.regions = colors(100), 
          at = seq(-.99, .99, .05))



#### Precipitation analysis vs IOD #############################################
# Prepare ONI record and classifiy years as IOD plus (P), IOD minus (M) or
# normal (N)
oni$TypeClass <- "N"
oni$TypeClass[grep("E", oni$Type)] <- "X"
oni$TypeClass[grep("L", oni$Type)] <- "X"
oni$TypeClass[grep("P", oni$IOD)] <- "IP"
oni$TypeClass[grep("M", oni$IOD)] <- "IM"

oni.reshape <- do.call("rbind", lapply(seq_len(nrow(oni)), function(i) {
  st <- as.Date(paste0(substr(oni[i, 3], 1, 4), "-07-01"))
  nd <- as.Date(paste0(substr(oni[i, 3], 6, 9), "-06-01"))
  data.frame(oni[i, c(1:3, ncol(oni))], 
             oni = as.numeric(oni[i, 4:(ncol(oni)-1)]), 
             month = seq(st, nd, "month"))
}))

# Adjust precipitation data records to ENSO cycle logic (i.e. start in July and
# end in June) and combine precipitation and ONI data set
precip.shift06m <- precip[7:(nrow(precip)-6), ]

precip.shift06m.oni <- 
  merge(precip.shift06m, oni.reshape, all.x = TRUE, by.x = "ts", by.y = "month")

precip.shift06m.oni.split <- 
  split(precip.shift06m.oni, precip.shift06m.oni$TypeClass)

precip.shift06m.oni.split.median <- 
  foreach(i = precip.shift06m.oni.split) %do% {
    sapply(1:12, function(j) {
      median(i$P_RT_NRT[seq(j, nrow(i), 12)], na.rm = TRUE)
    })
  }

# Produce publication quality graph for mean seasonal cycles during El Nino,
# La Nina and normal situations. To smooth the cycle while not reducing the
# rainfall amounts significantly, use a spline prediction.
if(length(precip.shift06m.oni.split) == 5){
  colors <- c("blue", "lightblue", "red", "bisque", "black")
} else {
  colors <- c("blue", "red", "black")  
}
xhres <- seq(1, 12, 0.01)
at <- seq(1, 1101, 100)
labels <- c(seq(7, 12, 1), seq(1, 6, 1))
plot.precip.shift06m.oni.split.median <- 
  lapply(seq(precip.shift06m.oni.split.median), function(x){
    sp <- predict(smooth.spline(
      precip.shift06m.oni.split.median[[x]], spar=0.01), xhres)
    xnew <- as.factor(sp[[1]])
    factor(xnew, levels = c(seq(7,12,0.01), seq(1,6,0.01)))
    xyplot(sp[[2]] ~ xnew, type = "l", lwd = "2", 
           ylim = c(0, 200), col = colors[x],
           scale=list(x=list(at = at, labels = labels)),
           xlab = "Month", ylab = "Precipitation")
  })
plot.precip.shift06m.oni.split.median.all <- 
  Reduce("outLayer", plot.precip.shift06m.oni.split.median)

# Non-spline version
# plot.precip.shift06m.oni.split.median <- 
#   lapply(seq(precip.shift06m.oni.split.median), function(x){
#     xyplot(precip.shift06m.oni.split.median[[x]] ~ 
#              factor(c(7:12, 1:6), levels = c(7:12, 1:6)), type = "l", 
#            ylim = c(0, 200), col = colors[x])
#   })
# 
# plot.precip.shift06m.oni.split.median.all <- 
#   Reduce("outLayer", plot.precip.shift06m.oni.split.median)
# plot.precip.shift06m.oni.split.median.all <- 
#   Reduce("outLayer", plot.precip.shift06m.oni.split.median)

# Just for background information: the individual seasons which form the above
# combined El Nino, La Nina and normal situation spline.
precip.shift06m.oni.split.season <- lapply(precip.shift06m.oni.split, function(x){
  split(x, x$Season)
})
xhres <- seq(1, 12, 0.01)
at <- seq(1, 1101, 100)
labels <- c(seq(7, 12, 1), seq(1, 6, 1))
plot.precip.shift06m.oni.split.season <- 
  lapply(seq(precip.shift06m.oni.split.season), function(x){
    lapply(seq(precip.shift06m.oni.split.season[[x]]), function(i) {
      sp <- with(precip.shift06m.oni.split.season[[x]][[i]][!is.na(precip.shift06m.oni.split.season[[x]][[i]]$P_RT_NRT),], smooth.spline(P_RT_NRT, spar=0.01))
      sp <- predict(sp, xhres)
      xnew <- as.factor(sp[[1]])
      factor(xnew, levels = c(seq(7,12,0.01), seq(1,6,0.01)))
      xyplot(sp[[2]] ~ xnew, type = "l", lwd = "2", 
             ylim = c(0, 700), col = colors[x],
             scale=list(x=list(at = at, labels = labels)),
             xlab = "Month", ylab = "Precipitation")
    })})
plot.precip.shift06m.oni.split.season.all <- 
  Reduce("outLayer", lapply(plot.precip.shift06m.oni.split.season, function(x){
    Reduce("outLayer", x)
  }))


# Alternative boxplot visualization
## bwplots
precip.shift06m.oni$month <- substr(precip.shift06m.oni$ts, 6, 7)
precip.shift06m.oni.melt <- melt(precip.shift06m.oni, 
                                 id.vars = c("month", "TypeClass"),
                                 measure.vars = "P_RT_NRT")
bw.plot <- ggplot(data = precip.shift06m.oni.melt, 
                  aes(y = value, x = month, fill = TypeClass, 
                      dodge = TypeClass))
bw.plot + geom_boxplot()

# Create publication quality figures of correlations between ONI and 
# precipitation
years <- unique(substr(precip$ts, 1, 4))
precip.shift06m.mat.ssn_kz03k01 <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.oni[grep(i, precip.shift06m.oni$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$ssn_kz03k01[j]
  }
}

oni.reshape.shift06m.mat <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.oni[grep(i, precip.shift06m.oni$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$oni[j]
  }
}

precip.shift06m.ssn_kz03k01.oni.reshape.shift06m.cor <- 
  cor(oni.reshape.shift06m.mat, precip.shift06m.mat.ssn_kz03k01, use = "complete.obs", 
      method = "kendall")
for (r in 1:12) {
  for (c in 1:12) {
    if (r > c) {
      precip.shift06m.ssn_kz03k01.oni.reshape.shift06m.cor[r,c] <- NA
    }
  }
}

colors <- colorRampPalette(brewer.pal(9, "RdBu"))
levelplot(t(precip.shift06m.ssn_kz03k01.oni.reshape.shift06m.cor), col.regions = colors(100), 
          at = seq(-.99, .99, .05))


# test <- subset(precip.shift06m.oni, ts >= "1997-01-01" & ts <= "2000-12-31")
# findMaxCCF(as.numeric(scale(test$ssn_kz03k02)), test$oni, na.action = na.exclude)
# 
# foreach(i = list(c("-03-", "-04-", "-05-"), 
#                  c("-06-", "-07-", "-08-"), 
#                  c("-09-", "-10-", "-11-"), 
#                  c("-12-", "-01-", "-02-"))) %do% {
#   sub <- subset(precip.shift06m.oni, ts >= "1997-01-01" & ts <= "2000-12-31")                 
#   index <- substr(sub$ts, 5, 8) %in% i
#   test2 <- sub[index, ]
#   findMaxCCF(as.numeric(scale(test2$ssn_kz03k02)), test2$oni, na.action = na.exclude)
# }
# 
# foreach(i = c("-04-", "-08-", "-12-")) %do% {
#   sub <- subset(precip.shift06m.oni, ts >= "1975-01-01" & ts <= "2013-12-31")                 
#   index <- substr(sub$ts, 5, 8) %in% i
#   test2 <- sub[index, ]
#   findMaxCCF(as.numeric(scale(test2$ssn_kz03k02)), test2$oni, na.action = na.exclude)
# }

# oni$IODAgg <- "N"
# oni$IODAgg[grep("P", oni$IOD)] <- "P"
# oni$IODAgg[grep("M", oni$IOD)] <- "M"
# 
# oni.reshape <- do.call("rbind", lapply(seq_len(nrow(oni)), function(i) {
#   st <- as.Date(paste0(substr(oni[i, 3], 1, 4), "-07-01"))
#   nd <- as.Date(paste0(substr(oni[i, 3], 6, 9), "-06-01"))
#   data.frame(oni[i, c(1:3, (ncol(oni)-1), ncol(oni))], 
#              oni = as.numeric(oni[i, 4:(ncol(oni)-2)]), 
#              month = seq(st, nd, "month"))
# }))
# 
# precip.shift06m.oni <- 
#   merge(precip.shift06m, oni.reshape, all.x = TRUE, by.x = "ts", by.y = "month")
# 
# precip.shift06m.oni.split <- 
#   split(precip.shift06m.oni, precip.shift06m.oni$IODAgg)
# 
# precip.shift06m.oni.split.median <- 
#   foreach(i = precip.shift06m.oni.split) %do% {
#     sapply(1:12, function(j) {
#       median(i$P_RT_NRT[seq(j, nrow(i), 12)], na.rm = TRUE)
#     })
#   }
# 
# plot(precip.shift06m.oni.split[[3]][1:12, "P_RT_NRT"], type = "l", ylim = c(0, 350), 
#      col = "blue")
# lines(precip.shift06m.oni.split[[3]][13:36, "P_RT_NRT"], col = "red")
# lines(precip.shift06m.oni.split[[3]][37:48, "P_RT_NRT"], col = "orange")
# lines(precip.shift06m.oni.split.median[[3]], col = "black")
# 
# vectorHarmonics(precip.shift06m.oni.split[[3]]$P_RT_NRT, st = c(1990, 1), 
#                 nd = c(2002, 12))

