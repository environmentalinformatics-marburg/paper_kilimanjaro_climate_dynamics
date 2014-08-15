seasonalMean <- function(data = precip,
                         st, nd,
                         st.shift = 7,
                         nd.shift = -6,
                         timespan = 12,
                         prm = "P_RT_NRT",
                         fun = "median"){
  
  if(st.shift == 0) {st.shift <- 1}
  
  if(fun == "mean"){
    foreach(st = st, nd = nd, .combine = "rbind") %do% {
      st.process <- grep(st, data[, 1])[st.shift]
      nd.process <- grep(nd, data[, 1])[length(grep(nd, data[, 1])) + nd.shift]
      data.process <- data[st.process:nd.process, ]
      df <- data.frame(season = paste(st, nd, sep = "-"), 
                       p_dyn = data.process.median <- 
                         sapply(1:12, function(j) {
                           mean(
                             data.process[, grep(prm, colnames(data.process))[1]][
                               seq(j, nrow(data.process), 12)], na.rm = TRUE)}))
      if(timespan != 12) {
        df <- rbind(df, df[c(1:(timespan - 12)),])
      }
      return(df)
    }
  }  else {
    foreach(st = st, nd = nd, .combine = "rbind") %do% {
      st.process <- grep(st, data[, 1])[st.shift]
      nd.process <- grep(nd, data[, 1])[length(grep(nd, data[, 1])) + nd.shift]
      data.process <- data[st.process:nd.process, ]
      df <- data.frame(season = paste(st, nd, sep = "-"), 
                       p_dyn = data.process.median <- 
                         sapply(1:12, function(j) {
                           median(
                             data.process[, grep(prm, colnames(data.process))[1]][
                               seq(j, nrow(data.process), 12)], na.rm = TRUE)}))
      if(timespan != 12) {
        df <- rbind(df, df[c(1:(timespan - 12)),])
      }
      return(df)
    }
  }
  
  
}
