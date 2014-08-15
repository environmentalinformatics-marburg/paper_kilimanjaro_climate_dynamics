#' Calculate and plot a correlation matrix of a time series data set
#' 
#' @description
#' The function calculates and plots a correlation matrix using two columns
#' of a data frame along the time line.
#' 
#' @param df a data frame
#' @param x.prm the name of the x values' column in df
#' @param y.prm the name of the y valus' volumn in df
#' @param t.prm the name of the time column in df
#' @param x.lable the x-axis lables of the plot
#' @param y.lable the y-axis lables of the plot
#' @param lable.nbrs the time numbering alon the axis
#' @param method the correlation test method
#' @param p.thv the p value to mark significant correlations in the plot
#' @param plot.filepath the path and filename of the output graphics file
#' @param plot2file controls if plot is written to file or standard out
#' @param cycle.window the window for the creation of the anomalies
#' 
#' @return a correlation matrix
#' 
#' @seealso
#' \code{\link{anomalize}}, \code{\link{denoise}}
#' 
#' @export visCorPlotTimeSeries
#' 
#' @examples 
#' none
visCorPlotTimeSeries <- function(df,
                                 x.prm,
                                 y.prm,
                                 t.prm,
                                 x.lable = substr(x.prm, 1, 3),
                                 y.lable = substr(y.prm, 1, 3),
                                 lable.nbrs = c(1:12),
                                 method = "kendall",
                                 p.thv = 0.05,
                                 plot.filepath = "corPlot.tif",
                                 plot2file = FALSE,
                                 rt = FALSE,
                                 rp = FALSE,
                                 x.text = NULL,
                                 y.text = NULL,
                                 labels.text = NULL,
                                 colors.text = NULL,
                                 ...){
  
  years <- unique(substr(df[, grep(t.prm, colnames(df))], 1, 4))
  
  x.data <- foreach(i = years, .combine = "rbind") %do% {
    dat <- df[grep(i, df[, grep(t.prm, colnames(df))]), ]
    foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
      dat[, grep(x.prm, colnames(dat))][j]
    }
  }
  colnames(x.data) <- paste0(x.lable, "-", sprintf("%02d", lable.nbrs))
  
  y.data <- foreach(i = years, .combine = "rbind") %do% {
    dat <- df[grep(i, df[, grep(t.prm, colnames(df))]), ]
    foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
      dat[, grep(y.prm, colnames(dat))[1]][j]
    }
  }
  colnames(y.data) <- paste0(y.lable, "-",sprintf("%02d", lable.nbrs))
  
  # Depricated
  cor.estimate <- cor(x.data, y.data, use = "pairwise.complete.obs", 
                      method = method)
  cor.estimate <- t(cor.estimate)
  
  cxy <- lapply(seq(ncol(x.data)), function(x) {
    vx <- x.data[,x]
    cy <- lapply(seq(ncol(y.data)), function(y) {
      vy <- y.data[,y]
      c <- cor.test(vx, vy, alternative = "two.sided", method = method)
      return(list(estimate = c$estimate,
                  pvalue = c$p.value))
    })
  })
  cor.estimate <- sapply(cxy, function(x) {
    sapply(x, function(y){y$estimate})})
  cor.pvalue <- sapply(cxy, function(x) {
    sapply(x, function(y){y$pvalue})})
  colnames(cor.estimate) <- colnames(x.data)
  rownames(cor.estimate) <- colnames(y.data)
  colnames(cor.pvalue) <- colnames(x.data)
  rownames(cor.pvalue) <- colnames(y.data)
  
  
  for (r in seq(nrow(cor.estimate))) {
    for (c in seq(ncol(cor.estimate))) {
      if (r < c) {
        cor.estimate[r,c] <- NA
      }
    }
  }
  
  cor.estimate.sig <- cor.estimate
  cor.estimate.sig[cor.pvalue > p.thv] <- NA
  print(max(cor.estimate.sig, na.rm = TRUE))
  print(min(cor.estimate.sig, na.rm = TRUE))
  
  colors <- colorRampPalette(brewer.pal(9, "RdBu"))
  plot.cor.estimate <- levelplot((cor.estimate), col.regions = colors(100), 
                                 at = seq(-.99, .99, .05), 
                                 scales = list(cex = 0.75, x = list(rot = 45)),
                                 xlab = "", ylab = "",
                                 colorkey = FALSE,
                                 panel = function(x,y,...){
                                   panel.levelplot(x,y,...)
                                   panel.text(x = x.text, y = y.text,
                                              labels = labels.text,
                                              col = colors.text,
                                              adj = c(0, 0))
                                 })
  plot.cor.estimate.sig <- levelplot((cor.estimate.sig), 
                                     col.regions = colors(100), 
                                     at = seq(-.99, .99, .05), 
                                     xlab = "", ylab = "", 
                                     border = "black")
  plot.cor <- doubleYScale(plot.cor.estimate, plot.cor.estimate.sig, 
                           add.axis = FALSE)
  
  if(plot2file == TRUE){
    tiff(filename = plot.filepath, 
         width = 2480, height = 1748 , res = 300, pointsize =  12)
    plot(plot.cor)
    dev.off()
  } else {
    plot(plot.cor)
  }
  if(rt){return(cor.estimate)}
  if(rp){return(plot.cor)}
}