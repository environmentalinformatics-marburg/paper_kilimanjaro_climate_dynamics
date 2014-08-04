seasonalMean <- function(st, nd){
  
  foreach(st = st, nd = nd, .combine = "rbind") %do% {
    st.process <- grep(st, precip[, 1])[7]
    nd.process <- grep(nd, precip[, 1])[length(grep(nd, precip[, 1]))-6]
    precip.process <- precip[st.process:nd.process, ]
    data.frame(season = paste(st, nd, sep = "-"), 
               p_dyn = precip.process.median <- 
                 sapply(1:12, function(j) {
                   median(precip.process$P_RT_NRT[seq(j, nrow(precip.process), 12)], na.rm = TRUE)
                 }))
  }
}