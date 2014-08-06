mergeAOIwTS <- function(aoi, 
                        ts,
                        ts.prm = "P_RT_NRT",
                        timespan = 12,
                        st.month = "-07-01",
                        nd.month = "-12-01",
                        rt = "median") {
  require(lubridate)
  aoi.reshape <- do.call("rbind", lapply(seq_len(nrow(aoi)), function(i) {
    st <- as.Date(paste0(substr(aoi[i, 3], 1, 4), st.month))
    nd <- st
    month(nd) <- month(nd) + (timespan - 1)
    aoi.prm = as.numeric(c(aoi[i, 4:(ncol(aoi)-1)],
                       aoi[i+1, 4:(ncol(aoi)-1)],
                       aoi[i+2, 4:(ncol(aoi)-1)]))
    aoi.prm = aoi.prm[1:timespan]
    data.frame(aoi[i, c(1:3, ncol(aoi))], 
               aoi = aoi.prm, 
               ts = as.character(seq(st, nd, "month")),
               stringsAsFactors=FALSE)
  }))
  aoi.reshape$SeasonCounter <- 
    paste0(rep((seq(1:(nrow(aoi.reshape)/timespan))+100), each = timespan), "-",
           aoi.reshape$ts)
  aoi.reshape$StartSeason <- substr(aoi.reshape$Season, 1, 4)
  ts$ts <- as.character(substr(ts$ts,1,10))
  ts.aoi <- 
    merge(ts, aoi.reshape, all.y = TRUE, by.x = "ts", by.y = "ts", 
          sort = FALSE)
  ts.aoi <- ts.aoi[with(ts.aoi, order(SeasonCounter)), ]
  
  ts.aoi <- ts.aoi[
    grep(TRUE, !is.na(
      ts.aoi[,grep(ts.prm, colnames(ts.aoi))[1]]))[1]:
      tail(grep(TRUE, !is.na(ts.aoi[,grep(ts.prm, colnames(ts.aoi))[1]])),1),]
  
  fst.complete.season <- 
    grep(as.character(as.numeric(ts.aoi$StartSeason[1]) + 1), 
         ts.aoi$StartSeason)[1]
  if(fst.complete.season >= timespan){
    fst.complete.season <- 1
  }
  
  lst.complete.season <- 
    tail(grep(as.character(as.numeric(tail(ts.aoi$StartSeason,1)) - 1),
              ts.aoi$StartSeason),1)
  if(nrow(ts.aoi) - lst.complete.season >= timespan){
    lst.complete.season <- nrow(ts.aoi)
  }
  
  ts.aoi <- ts.aoi[fst.complete.season:lst.complete.season,]

  ts.aoi.split <- 
    split(ts.aoi, ts.aoi$TypeClass)
  
  ts.aoi.split.median <- 
    foreach(i = ts.aoi.split) %do% {
      sapply(1:timespan, function(j) {
        median(i[, grep(ts.prm, colnames(i))[1]][seq(j, nrow(i), timespan)], na.rm = TRUE)
      })
    }
  
  ts.aoi.split.mean <- 
    foreach(i = ts.aoi.split) %do% {
      sapply(1:timespan, function(j) {
        mean(i[, grep(ts.prm, colnames(i))[1]][seq(j, nrow(i), timespan)], na.rm = TRUE)
      })
    }

  if(rt == "median"){
    return(ts.aoi.split.median)
  } else if(rt == "mean"){
    return(ts.aoi.split.mean)
  } else if(rt == "split") {
    return(ts.aoi.split)
  } else {
    return(ts.aoi)
  }
}
