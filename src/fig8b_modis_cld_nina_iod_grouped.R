library(reshape)
library(ggplot2)
library(RColorBrewer)
library(Rsenal)
library(grid)

# set working directory
setwd("C:/Users/IOtte/documents/Desktop/kilimanjaro/data/cloud_cover_monthly_mod06/cluster/") 

# read monthly precipitation data from Terra Modis
# combined with enso (La Nina) and IOD+
cld.nina.iod <- read.csv("clust_mean_series_terra_ensoiod.csv", 
                          header = TRUE, sep = ",")

### subset data
# overall cloud cover mean - 10 years
cld.nina.iod.all <- cld.nina.iod
cld.nina.iod.all$type <- "meanCld"


# all La Ninas - 5 years
cld.nina.iod.all.nina <- subset(cld.nina.iod, 
                           cld.nina.iod$ENSO == "ML" |
                           cld.nina.iod$ENSO == "SL" |
                           cld.nina.iod$ENSO == "WL")
cld.nina.iod.all.nina$type <- "all"
cld.nina.iod.alles <- rbind(cld.nina.iod.all, cld.nina.iod.all.nina)


# pure La Ninas - 4 years
cld.nina.iod.pure.enso <- subset(cld.nina.iod, 
                                 cld.nina.iod$ENSO == "ML" &
                                   cld.nina.iod$IOD != "P" &
                                   cld.nina.iod$IOD != "M" |  
                                   cld.nina.iod$ENSO == "SL" &
                                   cld.nina.iod$IOD != "P" &
                                   cld.nina.iod$IOD != "M" |  
                                   cld.nina.iod$ENSO == "WL" &
                                   cld.nina.iod$IOD != "P" &
                                   cld.nina.iod$IOD != "M")
cld.nina.iod.pure.enso$type <- "pureEnso"
cld.nina.iod.alles <- rbind(cld.nina.iod.alles, cld.nina.iod.pure.enso)


# La Nina IOD+ - 1 year
cld.nina.iod.enso.iod <- subset(cld.nina.iod, 
                                cld.nina.iod$ENSO == "ML" &
                                  cld.nina.iod$IOD == "P" |  
                                  cld.nina.iod$ENSO == "SL" &
                                  cld.nina.iod$IOD == "P" |
                                  cld.nina.iod$ENSO == "WL" &
                                  cld.nina.iod$IOD == "P")
cld.nina.iod.enso.iod$type <- "EnsoIOD"
cld.nina.iod.alles <- rbind(cld.nina.iod.alles, cld.nina.iod.enso.iod)


# purest IOD+ (incl. weak El Nino & La Nina) - 1 year
#cld.nina.iod.pure.iod <- subset(cld.nina.iod, 
#                                cld.nina.iod$IOD == "P" &
#                                cld.nina.iod$ENSO == "WL" |
#                                cld.nina.iod$IOD == "P" &
#                                cld.nina.iod$ENSO == "WE")
#cld.nina.iod.pure.iod$type <- "purestIOD"
#cld.nina.iod.alles <- rbind(cld.nina.iod.alles, cld.nina.iod.pure.iod)



## build factors for plotting boxplots
cld.nina.iod.alles$type <- factor(cld.nina.iod.alles$type, 
                                  levels = c("meanCld",
                                             "all", 
                                             "pureEnso",
                                             "EnsoIOD"))

# build factors for plotting boxplots
cld.nina.iod.alles$period <- factor(cld.nina.iod.alles$period, 
                                    levels = c("1", "2", "3", "4", "5", "6", 
                                               "7", "8", "9", "10", "11", "12",
                                               "13", "14", "15", "16", "17", 
                                               "18", "19", "20"))


## plot publication quality graphic of ENSO/IOD+ vs cloud clusters
### meanKIA/Moshi gemeinsam mit KIA und Moshi in einen plot darstellen
# melt cld.nina.iod.alles for facetting
cld.nina.iod.alles.mlt <- melt(cld.nina.iod.alles)
#copy
cld.nina.iod.alles.mlt.cp <- cld.nina.iod.alles.mlt
# rename type for facetting
levels(cld.nina.iod.alles.mlt.cp$variable)[levels(cld.nina.iod.alles.mlt.cp$variable) 
                                           == "clust1"] <- "cloud cluster a"
levels(cld.nina.iod.alles.mlt.cp$variable)[levels(cld.nina.iod.alles.mlt.cp$variable) 
                                           == "clust2"] <- "cloud cluster b"
levels(cld.nina.iod.alles.mlt.cp$variable)[levels(cld.nina.iod.alles.mlt.cp$variable) 
                                           == "clust3"] <- "cloud cluster c"

# subset (remove clust2, as the correlation is not significant)
cld.nina.iod.alles.mlt.cp.13 <- subset(cld.nina.iod.alles.mlt.cp, 
                                       cld.nina.iod.alles.mlt.cp$variable != 
                                         "cloud cluster b")

enso.iod <- subset(cld.nina.iod.alles.mlt.cp, cld.nina.iod.alles.mlt.cp$type == "EnsoIOD")
enso.iod$type <- factor(enso.iod$type, levels = c("meanCld", "all", "pureEnso", "EnsoIOD"))

cld.enso <- subset(cld.nina.iod.alles.mlt.cp, cld.nina.iod.alles.mlt.cp$type != "EnsoIOD")
cld.enso$type <- factor(cld.enso$type, levels = c("meanCld", "all", "pureEnso", "EnsoIOD"))


# dunkelgrau, blau, hellgrau, gr??n
cols <- c("#737373", brewer.pal(4, "Paired")[2], "#d9d9d9", brewer.pal(4, "Paired")[4])

# plotting
cld.nina.iod.boxplot.kia.moshi.mean.fct.sf <- ggplot(aes(x = period, y = value, 
                                                         fill = type), 
                                                     data = cld.enso) + 
  #geom_boxplot() +
  #facet_grid(variable ~.)
  geom_vline(xintercept = 6, color = "grey85") +
  geom_boxplot(position = "dodge", outlier.colour = "black", outlier.shape = 21) +
  geom_point(aes(x = period, y = value), data = enso.iod,
             position = "dodge", space = "free", size = 4, 
             color = brewer.pal(4, "Paired")[2], shape = 18) +
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
    legend.text = element_text(size = rel(.8)),
    legend.key.size = unit(.68, "cm"),
    legend.background = element_rect(color = "black", fill = "white"),
    #legend.position = c(.86,.908)) +       # legend position for aqua 
    legend.position = c(.65,.267)) +          # legend position for terra
  labs(fill = "") +
  facet_grid(variable ~ .) +
  scale_fill_manual(values = cols,
                    labels = c("all La Ninas (n=5)", 
                               "La Ninas w IOD+ (n=1)",
                                "mean cloud cover (n=10)",
                               "pure La Ninas (n=2)")) +
  scale_x_discrete(labels = c("07","08","09","10","11","12","01","02","03","04",
                              "05","06","07", "08", "09", "10", "11", "12","01", 
                              "02")) +
  xlab("Month") +
  ylab("Terra MODIS derived cloud cover [%]") 

# print cld.nina.iod.boxplot.kia.moshi.mean.fct.sf.png
png("out/cld.nina.iod.boxplot.kia.moshi.mean.fct.sf.aqua.png", width = 20, 
    height = 30, units = "cm", res = 300, pointsize = 15)
print(cld.nina.iod.boxplot.kia.moshi.mean.fct.sf)
dev.off()

