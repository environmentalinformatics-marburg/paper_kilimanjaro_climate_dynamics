visSeasonPlotByAOI <- function(data, colors, linetype = NULL,
                            individual = FALSE,
                            normal = NULL,
                            ymin = 0, ymax = 200,
                            timespan = 18,
                            xlab = "Month", ylab = "Precipitation",
                            vline.pos = NULL){
  
  xhres <- seq(1, timespan, 0.01)
  at <- seq(1, ((timespan-1)*100 + 1), 100)
  labels <- c(seq(7, 12, 1), rep(seq(1, 12, 1), ceiling((timespan-6)/12)))
  labels <- labels[1:timespan]
  if(is.null(linetype)){
    linetype = rep(1,length(data))
  }
  if(individual == TRUE){
    plot.data <- 
      lapply(seq(data), function(x){
        lapply(seq(data[[x]]), function(i) {
          sp <- with(data[[x]][[i]][!is.na(data[[x]][[i]]$P_RT_NRT),], 
                     smooth.spline(P_RT_NRT, spar=0.01))
          sp <- predict(sp, xhres)
          xnew <- as.factor(sp[[1]])
          factor(xnew, levels = c(seq(7,12,0.01), seq(1,6,0.01), seq(13,timespan,0.01)))
          xyplot(sp[[2]] ~ xnew, type = "l", lwd = "2", lty = linetype[x],
                 ylim = c(ymin, ymax), col = colors[x],
                 scale=list(x=list(at = at, labels = labels)),
                 xlab = xlab, ylab = ylab)
        })})
    plot.data.all <- 
      Reduce("outLayer", lapply(plot.data, function(x){
        Reduce("outLayer", x)
      }))    
  } else {
    plot.data <- 
      lapply(seq(data), function(x){
        sp <- predict(smooth.spline(
          data[[x]], spar=0.01), xhres)
        xnew <- as.factor(sp[[1]])
        factor(xnew, levels = c(seq(7,12,0.01), seq(1,6,0.01), seq(13,timespan,0.01)))
        xyplot(sp[[2]] ~ xnew, type = "l", lwd = "2", lty = linetype[x],
               ylim = c(ymin, ymax), col = colors[x],
               scale=list(x=list(at = at, labels = labels)),
               xlab = xlab, ylab = ylab, panel = function(x,y,...){
                 panel.xyplot(x,y,...)
                 panel.abline(v = vline.pos)
               })
      })
    plot.data.all <- 
      Reduce("outLayer", plot.data)
    plot.data.all <- plot.data.all
  }
  if(!is.null(normal)){
    plot.data.all <- 
      Reduce("outLayer", c(list(normal), plot.data))
  }
  return(plot.data.all)
}
