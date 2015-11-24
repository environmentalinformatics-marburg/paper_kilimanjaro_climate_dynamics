library(reshape)
library(ggplot2)
library(RColorBrewer)
library(Rsenal)

# set working directory
setwd("C:/Users/IOtte/documents/Desktop/kilimanjaro/data/metoffice/data") 

# read monthly precipitation data from the meteorological agency Tanzania
# combined with enso (La Nina) and iod
pre.nina.iod <- read.csv("metoffice_1973-2013_ensoiod_manipulated.csv", 
                          header = TRUE, sep = " ")

### subset data
# all precipitation data - 40 years
pre.nina.iod.all <- pre.nina.iod
pre.nina.iod.all$type <- "meanPre"


# all La Ninas - 15 years
pre.nina.iod.all.nina <- subset(pre.nina.iod, 
                                pre.nina.iod$ENSO == "ML" |
                                pre.nina.iod$ENSO == "SL" |
                                pre.nina.iod$ENSO == "WL")
pre.nina.iod.all.nina$type <- "all"
pre.nina.iod.alles <- rbind(pre.nina.iod.all, pre.nina.iod.all.nina)


# pure La Ninas - 11 years
pre.nina.iod.pure.enso <- subset(pre.nina.iod, 
                                 pre.nina.iod$ENSO == "ML" &
                                 pre.nina.iod$IOD != "P" &
                                 pre.nina.iod$IOD != "M" |  
                                 pre.nina.iod$ENSO == "SL" &
                                 pre.nina.iod$IOD != "P" &
                                 pre.nina.iod$IOD != "M" |  
                                 pre.nina.iod$ENSO == "WL" &
                                 pre.nina.iod$IOD != "P" &
                                 pre.nina.iod$IOD != "M")
pre.nina.iod.pure.enso$type <- "pureEnso"
pre.nina.iod.alles <- rbind(pre.nina.iod.alles, pre.nina.iod.pure.enso)


# La Nina IOD+ - 2 years
pre.nina.iod.enso.iod <- subset(pre.nina.iod, 
                                pre.nina.iod$ENSO == "ML" &
                                pre.nina.iod$IOD == "P" |  
                                pre.nina.iod$ENSO == "SL" &
                                pre.nina.iod$IOD == "P" |
                                pre.nina.iod$ENSO == "WL" &
                                pre.nina.iod$IOD == "P")
pre.nina.iod.enso.iod$type <- "EnsoIOD"
pre.nina.iod.alles <- rbind(pre.nina.iod.alles, pre.nina.iod.enso.iod)


# purest IOD+ (incl. weak El Nino) - 1 year
#pre.nina.iod.pure.iod <- subset(pre.nina.iod, 
#                                pre.nina.iod$IOD == "P" &
#                                pre.nina.iod$ENSO == "WL" |
#                                pre.nina.iod$IOD == "P" &
#                                pre.nina.iod$ENSO == "WE")
#pre.nina.iod.pure.iod$type <- "purestIOD"
#pre.nina.iod.alles <- rbind(pre.nina.iod.alles, pre.nina.iod.pure.iod)



## build factors for plotting boxplots
pre.nina.iod.alles$type <- factor(pre.nina.iod.alles$type, 
                                  levels = c("meanPre",
                                             "all", 
                                             "pureEnso",
                                             "EnsoIOD"))

# build factors for plotting boxplots
pre.nina.iod.alles$period <- factor(pre.nina.iod.alles$period, 
                                    levels = c("1", "2", "3", "4", "5", "6", 
                                               "7", "8", "9", "10", "11", "12",
                                               "13", "14", "15", "16", "17", 
                                               "18", "19", "20"))


## plot publication quality graphic of ENSO/IOD+ vs precipitation mean KIA/Moshi
### meanKIA/Moshi gemeinsam mit KIA und Moshi in einen plot darstellen
# melt pre.enso.iod.alles for facetting
pre.nina.iod.alles.mlt <- melt(pre.nina.iod.alles)
#copy
pre.nina.iod.alles.mlt.cp <- pre.nina.iod.alles.mlt
# rename type for facetting
levels(pre.nina.iod.alles.mlt.cp$variable)[levels(pre.nina.iod.alles.mlt.cp$variable) 
                                           == "P_KIA_NEW"] <- "KIA"
levels(pre.nina.iod.alles.mlt.cp$variable)[levels(pre.nina.iod.alles.mlt.cp$variable) 
                                           == "P_MOSHI_NEW"] <- "Moshi"
levels(pre.nina.iod.alles.mlt.cp$variable)[levels(pre.nina.iod.alles.mlt.cp$variable) 
                                           == "P_MEAN_NEW"] <- "mean KIA/Moshi"

cols <- c("#d9d9d9", "#737373", brewer.pal(4, "Paired")[c(4, 2)])


# plotting
pre.nina.iod.boxplot.kia.moshi.mean.fct.sf <- ggplot(aes(x = period, y = value, 
                                                  fill = type), 
                                              data = pre.nina.iod.alles.mlt.cp) + 
  #geom_boxplot() +
  #facet_grid(variable ~., scales = "free")
  geom_vline(xintercept = 6, color = "grey85") +
  geom_boxplot(position = "dodge", outlier.colour = "black", outlier.shape = 21) +
  facet_grid(variable ~ ., scales = "free_y") +
  scale_fill_manual(values = cols,
                    labels = c("mean rainfall (n=40)",
                               "all La Ninas (n=15)", 
                               "pure La Ninas (n=11)",
                               "La Ninas w IOD+ (n=2)")) +
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
    legend.position = c(.863,.596)) +
  labs(fill = "")

# print pre.enso.iod.boxplot.kia.moshi.mean.fct.sf.png
png("out/pre.nina.iod.boxplot.kia.moshi.mean.fct.sf.grpd.png", width = 20, 
    height = 30, units = "cm", res = 300, pointsize = 15)
print(pre.nina.iod.boxplot.kia.moshi.mean.fct.sf)
dev.off()

