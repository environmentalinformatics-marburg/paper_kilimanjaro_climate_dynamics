visSeasonPlotByAOI <- function(precip, colors, linetype = NULL,
                            individual = FALSE,
                            normal = NULL,
                            ymin = 0, ymax = 200,
                            xlab = "Month", ylab = "Precipitation"){
  
  xhres <- seq(1, 18, 0.01)
  at <- seq(1, 1701, 100)
  labels <- c(seq(7, 12, 1), seq(1, 12, 1))
  if(is.null(linetype)){
    linetype = rep(1,length(precip))
  }
  if(individual == TRUE){
    plot.precip <- 
      lapply(seq(precip), function(x){
        lapply(seq(precip[[x]]), function(i) {
          sp <- with(precip[[x]][[i]][!is.na(precip[[x]][[i]]$P_RT_NRT),], 
                     smooth.spline(P_RT_NRT, spar=0.01))
          sp <- predict(sp, xhres)
          xnew <- as.factor(sp[[1]])
          factor(xnew, levels = c(seq(7,12,0.01), seq(1,6,0.01), seq(13,18,0.01)))
          xyplot(sp[[2]] ~ xnew, type = "l", lwd = "2", lty = linetype[x],
                 ylim = c(ymin, ymax), col = colors[x],
                 scale=list(x=list(at = at, labels = labels)),
                 xlab = xlab, ylab = ylab)
        })})
    plot.precip.all <- 
      Reduce("outLayer", lapply(plot.precip, function(x){
        Reduce("outLayer", x)
      }))    
  } else {
    plot.precip <- 
      lapply(seq(precip), function(x){
        sp <- predict(smooth.spline(
          precip[[x]], spar=0.01), xhres)
        xnew <- as.factor(sp[[1]])
        factor(xnew, levels = c(seq(7,12,0.01), seq(1,6,0.01), seq(13,18,0.01)))
        xyplot(sp[[2]] ~ xnew, type = "l", lwd = "2", lty = linetype[x],
               ylim = c(ymin, ymax), col = colors[x],
               scale=list(x=list(at = at, labels = labels)),
               xlab = xlab, ylab = ylab)
      })
    plot.precip.all <- 
      Reduce("outLayer", plot.precip)
  }
  if(!is.null(normal)){
    plot.precip.all <- 
      Reduce("outLayer", c(list(normal), plot.precip))
  }
  return(plot.precip.all)
}
