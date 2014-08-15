longTermDynamicsPlot <- function(parameter = "temperature",
                                 p.prm = "ssn_kz03k02",
                                    printToFile = FALSE) {
  
  if(parameter == "temperature") {
    # Create publication quality figure of long term taitation trends
    greys <- brewer.pal(4, "Greys")
    colors <- c(greys[3], "red", "blue", greys[4])
    ylim = c(10,35)
    plot.ta.Ta_200 <- 
      xyplot(TEMP ~ Datetime, data = ta, origin = 0, type = "l",
             border = "transparent", asp = 0.25, col = colors[1],
             xlab = "", ylab = "Temperature (Deg C)", 
             lwd = 1.5, ylim = ylim, as.table = TRUE,
             scales = list(x = list(axs = "i")),
             xscale.components = xscale.components.subticks,
             yscale.components = yscale.components.subticks,
             panel = function(x, y, ...) {
               panel.xyplot(x, y, ...)
               panel.smoother(x, y, method = "lm", 
                              col = "black", 
                              col.se = "black",
                              alpha.se = 0.3, lty = 2, lwd = 2)
             })

    plot.ta.Ta_200.org <- 
      xyplot(TEMP ~ Datetime, data = ta.org, origin = 0, type = "l",
             border = "transparent", asp = 0.25, col = colors[4],
             xlab = "", ylab = "Temperature (Deg C)", 
             lwd = 1.5, ylim = ylim, as.table = TRUE,
             scales = list(x = list(axs = "i")),
             xscale.components = xscale.components.subticks,
             yscale.components = yscale.components.subticks,
             panel = function(x, y, ...) {
               panel.xyplot(x, y, ...)
               panel.smoother(x, y, method = "lm", 
                              col = "black", 
                              col.se = "black",
                              alpha.se = 0.3, lty = 2, lwd = 2)
             })
    
    plot.ta.Ta_200_Max <- 
      xyplot(MAX ~ Datetime, data = ta, origin = 0, type = "l",
             border = "transparent", asp = 0.25, col = colors[2],
             xlab = "", ylab = "Temperature ([Deg C)", 
             lwd = 1.5, ylim = ylim, as.table = TRUE,
             scales = list(x = list(axs = "i")),
             xscale.components = xscale.components.subticks,
             yscale.components = yscale.components.subticks,
             panel = function(x, y, ...) {
               #panel.xyplot(x, y, ...)
               panel.smoother(x, y, method = "lm", 
                              col = "red", 
                              col.se = "red",
                              alpha.se = 0.3, lty = 2, lwd = 2)
             })
    
    plot.ta.Ta_200_Min <- 
      xyplot(MIN ~ Datetime, data = ta, origin = 0, type = "l",
             border = "transparent", asp = 0.25, col = colors[2],
             xlab = "", ylab = "Temperature ([Deg C)", 
             lwd = 1.5, ylim = ylim, as.table = TRUE,
             scales = list(x = list(axs = "i")),
             xscale.components = xscale.components.subticks,
             yscale.components = yscale.components.subticks,
             panel = function(x, y, ...) {
               #panel.xyplot(x, y, ...)
               panel.smoother(x, y, method = "lm", 
                              col = "blue", 
                              col.se = "blue",
                              alpha.se = 0.3, lty = 2, lwd = 2)
             })
    
    plot.ta.Ta_200.all <- 
      Reduce("outLayer", list(plot.ta.Ta_200, plot.ta.Ta_200.org,
                              plot.ta.Ta_200_Max, plot.ta.Ta_200_Min))
    
    print(MannKendall(ta$TEMP))
    
    if(printToFile == TRUE){
      tiff(filename = paste0(graphicsPath, "plot.ta.Ta_200.all.tif"),
           width = 2480, height = 1748 , res = 300, pointsize =  12)
      plot(plot.ta.Ta_200.all)
      dev.off()
    } else {
      plot(plot.ta.Ta_200.all)  
    }
    
  } else if (parameter == "precipitation") {
    # Create publication quality figure of long term precipitation trends
    clr <- as.character(ifelse(precip[, grep(p.prm, colnames(precip))] > 0, "blue", "red"))
    clr.blue <- brewer.pal(9, "Blues")[2]
    clr.red <- brewer.pal(9, "Reds")[1]
    plot.precip.ssn_kz03k02 <- 
      xyplot(precip[, grep(p.prm, colnames(precip))] ~ precip$ts, origin = 0, type = "h",
             border = "transparent", col = clr, asp = 0.25,
             xlab = "", ylab = "Precipitation (mm)", 
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
    print(MannKendall(precip$ssn_kz03k02))
    
    if(printToFile == TRUE){
      tiff(filename = paste0(graphicsPath, "plot.precip.ssn_kz03k02.tif"),
           width = 2480, height = 1748 , res = 300, pointsize =  12)
      plot(plot.precip.ssn_kz03k02)
      dev.off()
    } else {
      plot(plot.precip.ssn_kz03k02)  
    }
    
  }  
}