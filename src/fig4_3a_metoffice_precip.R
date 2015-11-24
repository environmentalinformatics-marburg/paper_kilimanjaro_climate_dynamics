library(reshape)
library(ggplot2)
library(RColorBrewer)
library(Rsenal)

# set working directory
setwd("C:/Users/IOtte/documents/Desktop/kilimanjaro/data/metoffice/data/") 

# read monthly precipitation data from the meteorological agency Tanzania
precip <- read.csv("metoffice_1973-2013_copy.csv", header = TRUE, sep = ";")

# reformate date for analysis
precip$YEAR <- as.Date(precip$YEAR)

precip$mnth <- substr(precip$YEAR, 6, 7)

#####
# build yeraly sums
precip$ann <- substr(precip$YEAR, 1, 4)
precip.ann <- aggregate(precip[, 6], by = list(precip$ann), FUN = "sum")
colnames(precip.ann) <- c("Year", "rainfall_sum")

# plot yearly rainfall sums
precip.ann.pl <- ggplot(precip.ann, aes(x = Year, y = rainfall_sum)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", fill = "blue") +
  xlab("Year") +
  ylab("Precipitation") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.background = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA))

png("out/precip.ann.pl.png", width = 30, height = 22, units = "cm", 
    res = 300, pointsize = 15)
print(precip.ann.pl)
dev.off()

## cut data for analysis as the year starts in July and ends in June
## therefore cutting the first 6 months and the last 6 months

# cut first 6 months
precip.ct <- precip[-6, ]
precip.ct <- precip.ct[-5, ]
precip.ct <- precip.ct[-4, ]
precip.ct <- precip.ct[-3, ]
precip.ct <- precip.ct[-2, ]
precip.ct <- precip.ct[-1, ]

# cut last 6 months
precip.ct <- precip.ct[-486, ]
precip.ct <- precip.ct[-485, ]
precip.ct <- precip.ct[-484, ]
precip.ct <- precip.ct[-483, ]
precip.ct <- precip.ct[-482, ]
precip.ct <- precip.ct[-481, ]


## SE-Test
#subset1973-1983 january
#jan7383 <- subset(precip7383, precip7383$mnth == "01")
#jan7383mn <- mean(jan7383$P_MEAN_NEW)
#jan7383se <- se(jan7383$P_MEAN_NEW)
#jan7383sd <- sd(jan7383$P_MEAN_NEW)


### aggregate mean precipitation (KIA and Moshi) and calculate standard error
#1973-2013
precip.mnth <- aggregate(precip.ct[, 6], by = list(precip.ct$mnth), 
                         FUN = "mean")
colnames(precip.mnth) <- c("month", "mean7313")

precip.mnth.se <- aggregate(precip.ct[, 6], by = list(precip.ct$mnth), 
                            FUN = "se")
precip.mnth$se7313 <- precip.mnth.se[,2]

#1973-1983
precip7383 <- subset(precip.ct, precip.ct$YEAR < "1983-07-01")                                     
vrbl <- aggregate(precip7383[, 6], by = list(precip7383$mnth), FUN = "mean")
precip.mnth$mean7383 <- vrbl[,2]
vrbl <- aggregate(precip7383[, 6], by = list(precip7383$mnth), FUN = "se")
precip.mnth$se7383 <- vrbl[,2]

#1983-1993
precip8393 <- subset(precip.ct, precip.ct$YEAR > "1983-06-01" & precip.ct$YEAR < "1993-07-01")                                     
vrbl <- aggregate(precip8393[, 6], by = list(precip8393$mnth), FUN = "mean")
precip.mnth$mean8393 <- vrbl[,2]
vrbl <- aggregate(precip8393[, 6], by = list(precip8393$mnth), FUN = "se")
precip.mnth$se8393 <- vrbl[,2]

#1993-2003
precip9303 <- subset(precip.ct, precip.ct$YEAR > "1993-06-01" & precip.ct$YEAR < "2003-07-01")                                     
vrbl <- aggregate(precip9303[, 6], by = list(precip9303$mnth), FUN = "mean")
precip.mnth$mean9303 <- vrbl[,2]
vrbl <- aggregate(precip9303[, 6], by = list(precip9303$mnth), FUN = "se")
precip.mnth$se9303 <- vrbl[,2]

#2003-2013
precip0313 <- subset(precip.ct, precip.ct$YEAR > "2003-06-01")                                     
vrbl <- aggregate(precip0313[, 6], by = list(precip0313$mnth), FUN = "mean")
precip.mnth$mean0313 <- vrbl[,2]
vrbl <- aggregate(precip0313[, 6], by = list(precip0313$mnth), FUN = "se")
precip.mnth$se0313 <- vrbl[,2]

write.csv(precip.mnth, file = "precip.mnth.csv")

# reformate table
precip.mnth.mlt <- melt(precip.mnth)
colnames(precip.mnth.mlt) <- c("month", "mean_decade", "precipitation")

se <- subset(precip.mnth.mlt, substr(precip.mnth.mlt$mean_decade, 1, 2)
             == "se")

mean <- subset(precip.mnth.mlt, substr(precip.mnth.mlt$mean_decade, 1, 4)
               == "mean") 
mean$se <- se[, 3]


# plot barplot precipitation including errorbar
limits <- aes(ymax = mean$precipitation + mean$se, 
              ymin = mean$precipitation - mean$se)

precip.se.mean <- ggplot(mean, aes(x = month, y = precipitation, 
                                fill = mean_decade)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("black", brewer.pal(4, "Paired")[c(2,3,1,4)]),
                    labels = c("1973-2013", "1973-1983", "1983-1993", 
                               "1993-2003", "2003-2013")) +
  geom_errorbar(limits, position = "dodge", width = 0.9) +
  xlab("Month") +
  ylab("Precipitation") +
  xlim("07","08","09","10","11","12","01","02","03","04","05","06") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA),
    legend.position = c(.07,.87)) +
  labs(fill = "") 

# print seasonal_precip_barplot.png
png("out/seasonal_precipse_barplot.png", width = 30, height = 22, units = "cm", 
    res = 300, pointsize = 15)
print(precip.se.mean)
dev.off()


# plot stacked barplot
precip.mean.stck <- ggplot(mean, aes(x = month, y = precipitation, 
                                   fill = mean_decade)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("black", brewer.pal(4, "Paired")[c(2,3,1,4)]),
                    labels = c("1973-2013", "1973-1983", "1983-1993", 
                               "1993-2003", "2003-2013")) +
  geom_errorbar(limits) +
  xlab("Month") +
  ylab("Precipitation") +
  xlim("07","08","09","10","11","12","01","02","03","04","05","06") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA),
    legend.position = c(.07,.87)) +
  labs(fill = "") 

# print seasonal_precip_barplot.png
png("out/seasonal_precipstck_barplot.png", width = 30, height = 22, units = "cm", 
    res = 300, pointsize = 15)
print(precip.mean.stck)
dev.off()




### plot boxplot
precip.ct$season <- substr(precip.ct$YEAR, 1, 4)
precip.ct$season <- gsub("1992", "92", precip.ct$season)
write.csv(precip.ct, file = "precip.ct.csv")

precip.ct <- read.csv("precip.ct.csv", header = TRUE, sep = "\t")

colors <- c("#737373", brewer.pal(4, "Paired")[c(2,3,1,4)])

# test als subset
# subset 1973-1983
precip.ct$YEAR <- as.Date(precip.ct$YEAR)
#sub7383 <- subset(precip.ct, precip.ct$YEAR < "1983-07-01")

###


#precip.ct$mnth <- factor(precip.ct$mnth, levels = as.character(c(7:12, 1:6)))
precip.ct$type <- "all"
precip.ct.1973 <- subset(precip.ct, season == 1973)
precip.ct.1973$type <- "1973"
precip.ct.all <- rbind(precip.ct, precip.ct.1973)
precip.ct.1983 <- subset(precip.ct, season == 1983)
precip.ct.1983$type <- "1983"
precip.ct.all <- rbind(precip.ct.all, precip.ct.1983)
precip.ct.1993 <- subset(precip.ct, season == 1993)
precip.ct.1993$type <- "1993"
precip.ct.all <- rbind(precip.ct.all, precip.ct.1993)
precip.ct.2003 <- subset(precip.ct, season == 2003)
precip.ct.2003$type <- "2003"
precip.ct.all <- rbind(precip.ct.all, precip.ct.2003)

precip.ct.all$mnth <- factor(precip.ct.all$mnth, levels = as.character(c(7:12, 1:6)))
precip.ct.all$type <- factor(precip.ct.all$type, levels = c("all", "1973", "1983",
                                                            "1993", "2003"))

precip.boxplot <- ggplot(aes(x = mnth, y = P_MEAN_NEW, fill = type), 
                         data = precip.ct.all) + 
  scale_fill_manual(values = c("#737373", brewer.pal(4, "Paired")[c(2,3,1,4)]),
                      labels = c("1973-2013", "1973-1983", "1983-1993", 
                                 "1993-2003", "2003-2013")) +
  geom_boxplot(position = "dodge", outlier.colour = "black", outlier.shape = 21) +
  #stat_summary(fun.y = "mean", geom = "point", shape = 23, size = 3, position = "dodge") +
  xlab("Month") +
  ylab("Rainfall [mm]") +
  #theme(
  #  panel.grid.major = element_blank(),
  #  panel.grid.minor = element_blank(),
  #  panel.background = element_blank(),
  #  panel.border = element_rect(color = "grey", fill = NA),
  #  legend.position = c(.07,.87)) +
  #labs(fill = "")
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    strip.text = element_text(face = "bold", size = rel(1.5)),
    strip.background = element_rect(color = "black", fill = "white"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.background = element_rect(color = "black", fill = "white"),
    legend.position = c(.05,.885)) +
  labs(fill = "")


# print seasonal_precip_barplot.png
png("out/seasonal_precip_boxplot.png", width = 30, height = 22, units = "cm", 
    res = 300, pointsize = 15)
print(precip.boxplot)
dev.off()
