seasonalMean <- function(st, nd,
                         st.shift = 7,
                         nd.shift = -6){
  
  if(st.shift == 0) {st.shift <- 1}
  
  foreach(st = st, nd = nd, .combine = "rbind") %do% {
    st.process <- grep(st, precip[, 1])[st.shift]
    nd.process <- grep(nd, precip[, 1])[length(grep(nd, precip[, 1])) + nd.shift]
    precip.process <- precip[st.process:nd.process, ]
    df <- data.frame(season = paste(st, nd, sep = "-"), 
                     p_dyn = precip.process.median <- 
                       sapply(1:nd.pos, function(j) {
                         median(precip.process$P_RT_NRT[seq(j, nrow(precip.process), nd.pos)], na.rm = TRUE)
                       }))
    if(nd.shift == 0) {
      df <- rbind(df, df[c(1:6),])
    }
    return(df)
  }
}