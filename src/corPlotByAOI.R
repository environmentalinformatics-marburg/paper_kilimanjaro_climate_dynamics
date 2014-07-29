corPlotByAOI <- function(aoi.reshape.shift06m.mat,
                         precip.shift06m.mat.ssn_kz03k01){

  precip.shift06m.ssn_kz03k01.aoi.reshape.shift06m.cor <- 
    cor(aoi.reshape.shift06m.mat, precip.shift06m.mat.ssn_kz03k01, use = "complete.obs", 
        method = "kendall")
  
  print(max(precip.shift06m.ssn_kz03k01.aoi.reshape.shift06m.cor))
  print(min(precip.shift06m.ssn_kz03k01.aoi.reshape.shift06m.cor))
  for (r in 1:12) {
    for (c in 1:12) {
      if (r > c) {
        precip.shift06m.ssn_kz03k01.aoi.reshape.shift06m.cor[r,c] <- NA
      }
    }
  }
  
  colors <- colorRampPalette(brewer.pal(9, "RdBu"))
  levelplot(t(precip.shift06m.ssn_kz03k01.aoi.reshape.shift06m.cor), col.regions = colors(100), 
            at = seq(-.99, .99, .05), xlab = "", ylab = "")
}
