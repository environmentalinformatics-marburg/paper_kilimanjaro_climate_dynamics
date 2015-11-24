library(reshape)
library(ggplot2)
library(RColorBrewer)
library(Rsenal)
library(grid)


# set working directory
setwd("C:/Users/IOtte/documents/Desktop/kilimanjaro/data/cloud_cover_monthly_mod06/cluster/") 

# read monthly precipitation data from Terra Modis
# combined with enso (La Nina) and IOD+
cld.enso.iod <- read.csv("clust_mean_series_terra_ensoiod.csv", 
                          header = TRUE, sep = ",")

### subset data
# overall cloud cover mean
cld.enso.iod.all <- cld.enso.iod
cld.enso.iod.all$type <- "meanCld"


# all La Ninas - 3 years
cld.enso.iod.all.Nina <- subset(cld.enso.iod, 
                                cld.enso.iod$ENSO == "ME" |
                                cld.enso.iod$ENSO == "SE" |
                                cld.enso.iod$ENSO == "WE")
cld.enso.iod.all.Nina$type <- "all"
cld.enso.iod.alles <- rbind(cld.enso.iod.all, cld.enso.iod.all.Nina)


# pure La Ninas - 2 years
cld.enso.iod.pure.enso <- subset(cld.enso.iod, 
                                 cld.enso.iod$ENSO == "ME" &
                                   cld.enso.iod$IOD != "P" &
                                   cld.enso.iod$IOD != "M" |  
                                   cld.enso.iod$ENSO == "SE" &
                                   cld.enso.iod$IOD != "P" &
                                   cld.enso.iod$IOD != "M" |  
                                   cld.enso.iod$ENSO == "WE" &
                                   cld.enso.iod$IOD != "P" &
                                   cld.enso.iod$IOD != "M")
cld.enso.iod.pure.enso$type <- "pureEnso"
cld.enso.iod.alles <- rbind(cld.enso.iod.alles, cld.enso.iod.pure.enso)


# La Nina IOD+ - 1 year = purest IOD+
#cld.enso.iod.enso.iod <- subset(cld.enso.iod, 
#                                cld.enso.iod$ENSO == "ME" &
#                                  cld.enso.iod$IOD == "P" |  
#                                  cld.enso.iod$ENSO == "SE" &
#                                  cld.enso.iod$IOD == "P" |
#                                  cld.enso.iod$ENSO == "WE" &
#                                  cld.enso.iod$IOD == "P")
#cld.enso.iod.enso.iod$type <- "EnsoIOD"
#cld.enso.iod.alles <- rbind(cld.enso.iod.alles, cld.enso.iod.enso.iod)


# purest IOD+ (incl. weak El Nino & La Nina) - 1 year
cld.enso.iod.pure.iod <- subset(cld.enso.iod, 
                                cld.enso.iod$IOD == "P" &
                                cld.enso.iod$ENSO == "WL" |
                                cld.enso.iod$IOD == "P" &
                                cld.enso.iod$ENSO == "WE")
cld.enso.iod.pure.iod$type <- "purestIOD"
cld.enso.iod.alles <- rbind(cld.enso.iod.alles, cld.enso.iod.pure.iod)



## build factors for plotting boxplots
cld.enso.iod.alles$type <- factor(cld.enso.iod.alles$type, 
                                  levels = c("meanCld",
                                             "all", 
                                             "pureEnso",
                                             "purestIOD"))

# build factors for plotting boxplots
cld.enso.iod.alles$period <- factor(cld.enso.iod.alles$period, 
                                    levels = c("1", "2", "3", "4", "5", "6", 
                                               "7", "8", "9", "10", "11", "12",
                                               "13", "14", "15", "16", "17", 
                                               "18", "19", "20"))


## plot publication quality graphic of ENSO/IOD+ vs cloud clusters
### meanKIA/Moshi gemeinsam mit KIA und Moshi in einen plot darstellen
# melt cld.nina.iod.alles for facetting
cld.enso.iod.alles.mlt <- melt(cld.enso.iod.alles)
#copy
cld.enso.iod.alles.mlt.cp <- cld.enso.iod.alles.mlt
# rename type for facetting
levels(cld.enso.iod.alles.mlt.cp$variable)[levels(cld.enso.iod.alles.mlt.cp$variable) 
                                           == "clust1"] <- "cloud cluster a"
levels(cld.enso.iod.alles.mlt.cp$variable)[levels(cld.enso.iod.alles.mlt.cp$variable) 
                                           == "clust2"] <- "cloud cluster b"
levels(cld.enso.iod.alles.mlt.cp$variable)[levels(cld.enso.iod.alles.mlt.cp$variable) 
                                           == "clust3"] <- "cloud cluster c"

# subset (remove clust2, as the correlation is not significant)
cld.enso.iod.alles.mlt.cp.13 <- subset(cld.enso.iod.alles.mlt.cp, 
                                       cld.enso.iod.alles.mlt.cp$variable != 
                                         "cloud cluster b")

iod <- subset(cld.enso.iod.alles.mlt.cp, cld.enso.iod.alles.mlt.cp$type == "purestIOD")
cld.enso <- subset(cld.enso.iod.alles.mlt.cp, cld.enso.iod.alles.mlt.cp$type != "purestIOD")

cols <- c("#737373", "#d9d9d9", brewer.pal(4, "Paired")[4], "red")

# plotting
cld.enso.iod.boxplot.kia.moshi.mean.fct.sf <- ggplot(aes(x = period, y = value, 
                                                  fill = type), 
                                              data = cld.enso) + 
  #geom_boxplot() +
  #facet_grid(variable ~.)
  geom_vline(xintercept = 6, color = "grey85") +
  geom_boxplot(position = "dodge", outlier.colour = "black", outlier.shape = 21) +
  geom_point(aes(x = period, y = value), data = iod,
             position = "dodge", space = "free", size = 4, color = "red", shape = 18) +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    strip.text = element_text(face = "bold", size = rel(1.3)),
    strip.background = element_rect(color = "black", fill = "white"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.text = element_text(size = rel(.8)),
    legend.key.size = unit(.68, "cm"),
    legend.background = element_rect(color = "black", fill = "white"),
    #legend.position = c(.86,.938)) +       # legend position for aqua 
    legend.position = c(.65,.267)) +          # legend position for terra
  labs(fill = "") +
  facet_grid(variable ~ .) +
  scale_fill_manual(values = cols,
                    labels = c("all El Ninos (n=3)", "mean cloud cover (n=10)",
                               "pure El Ninos (n=2)","nearly pure IOD+ (n=1)")) +
  scale_x_discrete(labels = c("07","08","09","10","11","12","01","02","03","04",
                              "05","06","07", "08", "09", "10", "11", "12","01", 
                              "02")) +
  xlab("Month") +
  ylab("Aqua MODIS derived cloud cover [%]") 
  

# print cld.nina.iod.boxplot.kia.moshi.mean.fct.sf.png
png("out/cld.enso.iod.boxplot.kia.moshi.mean.fct.sf.grpd.terra.fake.png", width = 20, 
    height = 30, units = "cm", res = 300, pointsize = 15)
print(cld.enso.iod.boxplot.kia.moshi.mean.fct.sf)
dev.off()

