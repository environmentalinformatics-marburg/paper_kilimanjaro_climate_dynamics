library(reshape)
library(ggplot2)
library(RColorBrewer)
library(Rsenal)

# set working directory
setwd("C:/Users/IOtte/documents/Desktop/kilimanjaro/data/metoffice/data") 

# read monthly precipitation data from the meteorological agency Tanzania
# combined with enso and iod
pre.enso.iod <- read.csv("metoffice_1973-2013_ensoiod_manipulated.csv", 
                         header = TRUE, sep = " ")

### subset data
# all precipitation data - 40 years
pre.enso.iod.all <- pre.enso.iod
pre.enso.iod.all$type <- "meanPre"


# all El Nino - 12 years
pre.enso.iod.all.nino <- subset(pre.enso.iod, 
                           pre.enso.iod$ENSO == "ME" |
                           pre.enso.iod$ENSO == "SE" |
                           pre.enso.iod$ENSO == "WE")
pre.enso.iod.all.nino$type <- "all"
pre.enso.iod.alles <- rbind(pre.enso.iod.all, pre.enso.iod.all.nino)


# pure El Nino - 7 years
pre.enso.iod.pure.enso <- subset(pre.enso.iod, 
                                 pre.enso.iod$ENSO == "ME" &
                                 pre.enso.iod$IOD != "P" &
                                 pre.enso.iod$IOD != "M" |  
                                 pre.enso.iod$ENSO == "SE" &
                                 pre.enso.iod$IOD != "P" &
                                 pre.enso.iod$IOD != "M" |  
                                 pre.enso.iod$ENSO == "WE" &
                                 pre.enso.iod$IOD != "P" &
                                 pre.enso.iod$IOD != "M")
pre.enso.iod.pure.enso$type <- "pureEnso"
pre.enso.iod.alles <- rbind(pre.enso.iod.alles, pre.enso.iod.pure.enso)


# pure m/s EL Nino - 5 years
#pre.enso.iod.pure.enso.ms <- subset(pre.enso.iod, 
#                                    pre.enso.iod$ENSO == "ME" &
#                                    pre.enso.iod$IOD != "P" |  
#                                    pre.enso.iod$ENSO == "SE" &
#                                    pre.enso.iod$IOD != "P")
#pre.enso.iod.pure.enso.ms$type <- "pureEnsoMS"
#pre.enso.iod.alles <- rbind(pre.enso.iod.alles, pre.enso.iod.pure.enso.ms)


# El Nino IOD+ - 5 years
pre.enso.iod.enso.iod <- subset(pre.enso.iod, 
                                pre.enso.iod$ENSO == "ME" &
                                pre.enso.iod$IOD == "P" |  
                                pre.enso.iod$ENSO == "SE" &
                                pre.enso.iod$IOD == "P" |
                                pre.enso.iod$ENSO == "WE" &
                                pre.enso.iod$IOD == "P")
pre.enso.iod.enso.iod$type <- "EnsoIOD"
pre.enso.iod.alles <- rbind(pre.enso.iod.alles, pre.enso.iod.enso.iod)


# m/s El Nino IOD+ - 3 years
#pre.enso.iod.enso.iod.ms <- subset(pre.enso.iod, 
#                                   pre.enso.iod$ENSO == "ME" &
#                                   pre.enso.iod$IOD == "P" |
#                                   pre.enso.iod$ENSO == "SE" &
#                                   pre.enso.iod$IOD == "P")
#pre.enso.iod.enso.iod.ms$type <- "EnsoMSIOD"
#pre.enso.iod.alles <- rbind(pre.enso.iod.alles, pre.enso.iod.enso.iod.ms)


# purest IOD+ (incl. weak El Nino & weak La Nina) - 3 years
pre.enso.iod.pure.iod <- subset(pre.enso.iod, 
                                pre.enso.iod$IOD == "P" &
                                pre.enso.iod$ENSO == "WE"|
                                pre.enso.iod$IOD == "P" &
                                pre.enso.iod$ENSO == "WL")
pre.enso.iod.pure.iod$type <- "purestIOD"
pre.enso.iod.alles <- rbind(pre.enso.iod.alles, pre.enso.iod.pure.iod)



## build factors for plotting boxplots
pre.enso.iod.alles$type <- factor(pre.enso.iod.alles$type, 
                                  levels = c("meanPre",
                                             "all", 
                                             "pureEnso",
                                             "EnsoIOD",
                                             "purestIOD"))

# build factors for plotting boxplots
pre.enso.iod.alles$period <- factor(pre.enso.iod.alles$period, 
                                    levels = c("1", "2", "3",
                                               "4", "5", "6",
                                               "7", "8", "9", 
                                               "10", "11", "12",
                                               "13", "14", "15",
                                               "16", "17", "18", "19", "20"))



## plot publication quality graphic of ENSO/IOD+ vs precipitation mean KIA/Moshi
pre.enso.iod.boxplot.mean.kia.moshi <- ggplot(aes(x = period, y = P_MEAN_NEW, 
                                                  fill = type), 
                                              data = pre.enso.iod.alles) + 
  geom_boxplot(position = "dodge") +
  # purest IOD+ nicht als boxplot sondern nur 2 punkte
  #geom_point(aes(x = period, y = P_MEAN_NEW), data = pre.enso.iod.pure.iod,
  #           position = "dodge", size = 6, color = "red") +
  geom_vline(xintercept = 6) +
  scale_fill_manual(values = c("black", brewer.pal(4, "Paired")[c(2,1,4,3)], 
                               "red"),
                    labels = c("all El Ninos", 
                               "El Ninos w IOD+",
                               "m/s El Ninos w IOD+",
                               "pure El Ninos",
                               "pure m/s El Ninos",
                               "purest IOD+")) +
  scale_x_discrete(labels = c("07","08","09","10","11","12","01","02","03","04",
                              "05","06","07", "08", "09", "10", "11", "12","01", 
                              "02")) +
  #scale_y_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300, 350)) +
  xlab("Month") +
  ylab("Precipitation") +
  annotate("text", x = 2.1, y = 380, label = "c) mean KIA/Moshi") +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.background = element_rect(color = "black", fill = "white"),
    legend.position = c(.923,.85)) +
  labs(fill = "")

# print spre.enso.iod.boxplot.mean.kia.moshi.png
png("out/pre.enso.iod.boxplot.mean.kia.moshi.png", width = 30, height = 20, 
    units = "cm", res = 300, pointsize = 15)
print(pre.enso.iod.boxplot.mean.kia.moshi)
dev.off()



### meanKIA/Moshi gemeinsam mit KIA und Moshi in einen plot darstellen
# melt pre.enso.iod.alles for facetting
pre.enso.iod.alles.mlt <- melt(pre.enso.iod.alles)
#copy
pre.enso.iod.alles.mlt.cp <- pre.enso.iod.alles.mlt
# rename type for facetting
levels(pre.enso.iod.alles.mlt.cp$variable)[levels(pre.enso.iod.alles.mlt.cp$variable) 
                                           == "P_KIA_NEW"] <- "KIA"
levels(pre.enso.iod.alles.mlt.cp$variable)[levels(pre.enso.iod.alles.mlt.cp$variable) 
                                           == "P_MOSHI_NEW"] <- "Moshi"
levels(pre.enso.iod.alles.mlt.cp$variable)[levels(pre.enso.iod.alles.mlt.cp$variable) 
                                           == "P_MEAN_NEW"] <- "mean KIA/Moshi"


cols <- c("#d9d9d9", "#737373", brewer.pal(4, "Paired")[c(4, 2)], "red")


# plotting
pre.enso.iod.boxplot.kia.moshi.mean.fct.sf <- ggplot(aes(x = period, y = value, 
                                                  fill = type), 
                                              data = pre.enso.iod.alles.mlt.cp) + 
  #geom_boxplot() +
  #facet_grid(variable ~., scales = "free")
  geom_vline(xintercept = 6, color = "grey85") +
  geom_boxplot(position = "dodge", outlier.colour = "black", outlier.shape = 21) +
  # purest IOD+ nicht als boxplot sondern nur 2 punkte
  #geom_point(aes(x = period, y = P_MEAN_NEW), data = pre.enso.iod.pure.iod,
  #           position = "dodge", size = 6, color = "red") +
  facet_grid(variable ~ ., scales = "free_y") +
  scale_fill_manual(values = cols,
                    labels = c("mean rainfall (n=40)",
                               "all El Ninos (n=12)", 
                               "pure El Ninos (n=7)",
                               "El Ninos w IOD+ (n=5)",
                               "nearly pure IOD+ (n=3)")) +
  scale_x_discrete(labels = c("07","08","09","10","11","12","01","02","03","04",
                              "05","06","07", "08", "09", "10", "11", "12","01", 
                              "02")) +
  xlab("Month") +
  ylab("Rainfall [mm]") +
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
    legend.position = c(.864,.582)) +
  labs(fill = "")

# print pre.enso.iod.boxplot.kia.moshi.mean.fct.sf.png
png("out/pre.enso.iod.boxplot.kia.moshi.mean.fct.sf.grpd.png", width = 20, 
    height = 30, units = "cm", res = 300, pointsize = 15)
print(pre.enso.iod.boxplot.kia.moshi.mean.fct.sf)
dev.off()

#####test#####
## try to color the outliers seperatly
#library(plyr)
#plot_Data <- ddply(pre.enso.iod.alles.mlt.cp, .(period, type), mutate, 
#                   Q1=quantile(value, 1/4), Q3=quantile(value, 3/4), IQR=Q3-Q1, 
#                   upper.limit=Q3+1.5*IQR, lower.limit=Q1-1.5*IQR)

#ggplot() +
#  geom_boxplot(data=plot_Data, aes(x=factor(period), y=value, col=factor(type))) + 
#  geom_point(data=plot_Data[plot_Data$value > plot_Data$upper.limit | 
#                              plot_Data$value < plot_Data$lower.limit,], 
#             aes(x=factor(period), y=value, col=factor(type)))



### mean and se per subset f??r KIA
#precip.mnth <- aggregate(precip.ct[, 6], by = list(precip.ct$mnth), 
#                         FUN = "mean")
#colnames(precip.mnth) <- c("month", "mean7313")

#precip.mnth.se <- aggregate(precip.ct[, 6], by = list(precip.ct$mnth), 
#                            FUN = "se")
#precip.mnth$se7313 <- precip.mnth.se[,2]

# all El Nino
#enso.all <- aggregate(pre.enso.iod.all[, 2], 
#                      by = list(substr(pre.enso.iod.all$YEAR, 6, 7)),
#                                                       FUN = "mean")
#colnames(enso.all) <- c("month", "precip_enso_all")
#enso.all.se <- aggregate(pre.enso.iod.all[, 2], 
#                         by = list(substr(pre.enso.iod.all$YEAR, 6, 7)),
#                         FUN = "se")
#enso.all$se_precip_enso_all <- enso.all.se[, 2]

# pure El Nino
# pure m/s EL Nino
# El Nino IOD+
# m/s El Nino IOD+
# purest IOD+ (incl. weak El Nino)
