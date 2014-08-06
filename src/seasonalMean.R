seasonalMean <- function(st, nd,
                         st.shift = 7,
                         nd.shift = -6,
                         timespan = 12,
                         prm = "P_RT_NRT",
                         fun = "median"){
  
  if(st.shift == 0) {st.shift <- 1}
  
  if(fun == "mean"){
    foreach(st = st, nd = nd, .combine = "rbind") %do% {
      st.process <- grep(st, precip[, 1])[st.shift]
      nd.process <- grep(nd, precip[, 1])[length(grep(nd, precip[, 1])) + nd.shift]
      precip.process <- precip[st.process:nd.process, ]
      df <- data.frame(season = paste(st, nd, sep = "-"), 
                       p_dyn = precip.process.median <- 
                         sapply(1:12, function(j) {
                           mean(
                             precip.process[, grep(prm, colnames(precip.process))[1]][
                               seq(j, nrow(precip.process), 12)], na.rm = TRUE)}))
      if(timespan != 12) {
        df <- rbind(df, df[c(1:(timespan - 12)),])
      }
      return(df)
    }
  }  else {
    foreach(st = st, nd = nd, .combine = "rbind") %do% {
      st.process <- grep(st, precip[, 1])[st.shift]
      nd.process <- grep(nd, precip[, 1])[length(grep(nd, precip[, 1])) + nd.shift]
      precip.process <- precip[st.process:nd.process, ]
      df <- data.frame(season = paste(st, nd, sep = "-"), 
                       p_dyn = precip.process.median <- 
                         sapply(1:12, function(j) {
                           median(
                             precip.process[, grep(prm, colnames(precip.process))[1]][
                               seq(j, nrow(precip.process), 12)], na.rm = TRUE)}))
      if(timespan != 12) {
        df <- rbind(df, df[c(1:(timespan - 12)),])
      }
      return(df)
    }
  }
  
  
}
