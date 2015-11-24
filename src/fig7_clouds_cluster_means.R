library(TSA)
library(latticeExtra)
library(raster)
library(TSclust)
library(wesanderson)
library(fpc)
library(cluster)
library(remote)
library(gridBase)

setwd("C:/Users/IOtte/documents/Desktop/kilimanjaro/data/cloud_cover_monthly_myd06/cluster/") 


#source("cellHarmonics.R")

clouds <- stack("stck_aq_mnth_0313.tif")

### cluster spatially deseasoned
clouds_dsn <- deseason(clouds, 12)

clouds_dsn_df <- as.data.frame(clouds_dsn)
clusters_dsn <- clouds_dsn[[1]]

#clust_clouds_dsn_df <- clara(clouds_dsn_df, 4, metric = "euclidian", samples = 1000) 
clust_clouds_dsn_df2 <- kmeans(clouds_dsn_df, 3, nstart = 4, iter.max = 100, algorithm = "Ll") ### kmeans, with 4 clusters

clusters_dsn[] <- clust_clouds_dsn_df2$cluster

plot(clusters_dsn)
clusplot(clouds_dsn_df, clust_clouds_dsn_df2$cluster)
d <- dist(clouds_dsn_df)
cluster.stats(d, clustering = clust_clouds_dsn_df2$cluster)

### mean frequency per cluster
meanFreq <- function(x) {
  ind <- which(clusters_dsn[] == x)
  as.numeric(colMeans(clouds[][ind, ], na.rm = TRUE))
}

df_clust <- data.frame("clust1" = meanFreq(1),
                       "clust2" = meanFreq(2),
                       "clust3" = meanFreq(3))

write.csv(df_clust, "clust_mean_series_aqua_controll.csv", row.names = FALSE)



# plotting
clrs_clust <- colorRampPalette(brewer.pal(5, "PuBu")[c(5,4,2)])
#clrs_rsq <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")))
ki.dem.utm <- raster("DEM_UTM37S_WGS84_30m_Hemp.tif")

panel.smoothconts <- function(x, y, z, at, 
                              contours = T, 
                              zlevs.conts = seq(500, 6000, 500),
                              ...)
{
  stopifnot(require("gridBase"))
  z <- matrix(z,
              nrow = length(unique(x)),
              ncol = length(unique(y)))
  rownames(z) <- unique(coordinates(ki.dem.utm)[, 1])
  colnames(z) <- unique(coordinates(ki.dem.utm)[, 2])
  
  if (!is.double(z)) storage.mode(z) <- "double"
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  if (panel.number() > 1) par(new = TRUE)
  par(fig = gridFIG(), omi = c(0, 0, 0, 0), mai = c(0, 0, 0, 0))
  cpl <- current.panel.limits(unit = "native")
  plot.window(xlim = cpl$xlim, ylim = cpl$ylim,
              log = "", xaxs = "i", yaxs = "i")
  # paint the color contour regions
  
  if (isTRUE(contours)) 
    contour(as.double(do.breaks(range(as.numeric(rownames(z))), nrow(z) - 1)),
            as.double(do.breaks(range(as.numeric(colnames(z))), ncol(z) - 1)),
            z, levels = as.double(zlevs.conts), 
            add = T, cex = 1.8,
            axes = FALSE,
            col = "grey60", # color of the lines
            drawlabels = TRUE,
            labcex = 1 # add labels or not
    )
}

ki.dem.flipped <- flip(ki.dem.utm, "y")
x <- coordinates(ki.dem.flipped)[, 1]
y <- coordinates(ki.dem.flipped)[, 2]
z <- ki.dem.flipped[]

cp <- levelplot(z ~ x * y, alpha = 0.7,
                at = seq(500, 6000, 500), colorkey = FALSE,
                panel = function(...) {
                  panel.smoothconts(contours = TRUE, ...)
                })

panel_txt_clust <- "d) Kmeans clusters"

clust_p <- spplot(clusters_dsn, col.regions = clrs_clust(1000),
                  colorkey = FALSE, scales = list(draw = TRUE, 
                                                  y = list(rot = 90),
                                                  alternating = 3),
                  panel = function(x, y, ...) {
                    panel.levelplot(x, y, ...)
                    #                     panel.text(x = 280000, y = 9627500, 
                    #                                labels = panel_txt_clust,
                    #                                adj = c(0, 0.5),
                    #                                col = "black")
                  }) + 
  as.layer(cp, x.same = TRUE, y.same = TRUE)

png("clouds_clusters_try.png", width = 25, 
    height = 20, units = "cm", res = 300)
plot.new()
print(clust_p)
dev.off()
