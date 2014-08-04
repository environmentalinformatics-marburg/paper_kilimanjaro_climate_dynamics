readData <- function(parameter = "temperature") {
  
  if(parameter == "temperature") {
    taKIA <- 
      read.csv(
        paste0(dataPath, 
               "temperature/noaa_gsod_kia_1973_2013_reformat_sd_li_jl_ssa.csv"), 
        stringsAsFactors = FALSE)
    taKIA[, 1] <- as.Date(taKIA[, 1])
    taKIA$YearMonth <- substr(taKIA$Datetime,1,7)
    taKIAMean <- aggregate(taKIA$TEMP, by = list(taKIA$YearMonth), FUN = "mean")
    colnames(taKIAMean) <- c("ts","Ta_200")
    taKIAMean$ts <- as.Date(paste0(taKIAMean$ts, "-01"))
    taKIAMean$Ta_200_Max <- aggregate(
      taKIA$MAX, by = list(taKIA$YearMonth), FUN = "mean")[[2]]
    taKIAMean$Ta_200_Min <- aggregate(
      taKIA$MIN, by = list(taKIA$YearMonth), FUN = "mean")[[2]]
    return(list(KIA = taKIA, KIAMean = taKIAMean))
  
  } else if (parameter == "precipitation") {
    precipKIA <- 
      read.csv(
        paste0(dataPath, 
               "kilimanjaro_gsod_dynamics/gsod_precip/kia_prcp_1975_2014_mnthly.csv"), 
        stringsAsFactors = FALSE)
    precipKIA[, 1] <- as.Date(precipKIA[, 1])
    st <- precipKIA[1, 1]
    nd <- precipKIA[nrow(precipKIA), 1]
    ts <- seq(st, nd, "month")
    precipKIA <- merge(data.frame(ts), precipKIA, by = 1, all.x = TRUE)
    
    precipNRB <- 
      read.csv(
        paste0(dataPath, 
               "kilimanjaro_gsod_dynamics/gsod_precip/nairobi_prcp_1975_2014_mnthly.csv"), 
        stringsAsFactors = FALSE)
    precipNRB[, 1] <- as.Date(precipNRB[, 1])
    st <- precipNRB[1, 1]
    nd <- precipNRB[nrow(precipNRB), 1]
    ts <- seq(st, nd, "month")
    precipNRB <- merge(data.frame(ts), precipNRB, by = 1, all.x = TRUE)
    return(list(KIA = precipKIA, NRB = precipNRB))
  } else if(parameter == "cloudEOT") {
    cloudEOT <- read.csv(paste0(dataPath, "cloud_cover_monthly_myd06/eot/eot_series.csv"), 
                         header = TRUE, sep = ",")
    cloudEOT$Years <- rep(2003:2013, each=12)
    cloudEOT$Month <- rep(1:12)
    cloudEOT$ts <- as.Date(paste0(cloudEOT$Year,"-",cloudEOT$Month,"-1"))
    cloudEOT <- cloudEOT[,c(6,1:3)]
    
    # Read deseasoned cloud EOT
    cloudEOTssn <- read.csv(paste0(dataPath, "cloud_cover_monthly_myd06/eot/eot_series_ssn.csv"), 
                            header = TRUE, sep = ",")
    cloudEOTssn$Years <- rep(2003:2013, each=12)
    cloudEOTssn$Month <- rep(1:12)
    cloudEOTssn$ts <- as.Date(paste0(cloudEOTssn$Year,"-",cloudEOTssn$Month,"-1"))
    cloudEOTssn <- cloudEOTssn[,c(8,2:4)]
    return(list(EOT = cloudEOT, EOTssn = cloudEOTssn))
  } else if(parameter == "aoi") {
    # Read ENSO ONI index
    oni <- read.csv(paste0(dataPath, "aoi/enso_and_iod.csv"),
                    skip = 1, header = TRUE)
    oni$Season <- paste0(oni$Season, oni$X, oni$X.1)
    oni <- oni[, -grep("X", names(oni))]
    
    # Read ENSO MEI index
    mei <- read.csv(paste0(dataPath, "aoi/mei.txt"), skip = 9, header = TRUE, 
                    sep = "\t", nrows = 65)
    mei <- mei[,c(1,8:13,2:7)]
    mei <- do.call(rbind, lapply(seq(nrow(mei)-1), function(x){
      modrow <- mei[x,]
      modrow[,8:13] <- mei[x+1,8:13]
      return(modrow)
    }))
    mei$Season <- paste0(mei$YEAR,"-",mei$YEAR+1)
    mei <- cbind(oni[,1:3],mei[,2:13])
    
    # Read DMI index
    dmi <- read.csv(paste0(dataPath, "aoi/dmi_in.txt"), header = FALSE, sep = ":")
    colnames(dmi) <- c("Date", "DMI")
    dmi$YearMonth <- substr(dmi$Date,1,7)
    dmi <- aggregate(dmi$DMI, by=list(dmi$YearMonth), FUN = "mean")
    colnames(dmi) <- c("YearMonth", "DMI")
    dmi <- dmi[c((grep("1981",dmi$YearMonth)[-1]+1):(grep("2014-07",
                                                          dmi$YearMonth)[1]-1)),]
    years <- unique(substr(dmi$YearMonth, 1, 4))
    dmi <- as.data.frame(foreach(i = years, .combine = "rbind") %do% {
      dat <- dmi[grep(i, dmi$YearMonth), ]
      foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
        dat$DMI[j]
      }
    })
    dmi$Year <- years
    dmi <- dmi[,c(13,7:12,1:6)]
    dmi <- do.call(rbind, lapply(seq(nrow(dmi)-1), function(x){
      modrow <- dmi[x,]
      modrow[,8:13] <- dmi[x+1,8:13]
      return(modrow)
    }))
    dmi$Season <- paste0(as.numeric(dmi$Year),"-",as.numeric(dmi$Year)+1)
    fr <- grep(dmi$Season[1], oni$Season)
    lr <- grep(dmi$Season[nrow(dmi)], oni$Season)
    dmi <- cbind(oni[fr:lr,1:3],dmi[,2:13])
    colnames(dmi)[4:15] <- colnames(oni)[4:15]
    
    # Read DMI-sstHAD index
    dmiHAD <- read.csv(paste0(dataPath, "aoi/dmi_hadisst.txt"), 
                       header = FALSE, sep = " ")
    colnames(dmiHAD) <- c("Year", "Month", "DMIHAT")
    dmiHAD$YearMonth <- paste0(dmiHAD$Year,"-",sprintf("%02d", dmiHAD$Month))
    dmiHAD <- dmiHAD[,c(4,3)]
    years <- unique(substr(dmiHAD$YearMonth, 1, 4))
    dmiHAD <- as.data.frame(foreach(i = years, .combine = "rbind") %do% {
      dat <- dmiHAD[grep(i, dmiHAD$YearMonth), ]
      foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
        dat$DMIHAT[j]
      }
    })
    dmiHAD$Year <- years
    dmiHAD <- dmiHAD[,c(13,7:12,1:6)]
    dmiHAD <- do.call(rbind, lapply(seq(nrow(dmiHAD)-1), function(x){
      modrow <- dmiHAD[x,]
      modrow[,8:13] <- dmiHAD[x+1,8:13]
      return(modrow)
    }))
    dmiHAD$Season <- paste0(as.numeric(dmiHAD$Year),"-",as.numeric(dmiHAD$Year)+1)
    fr <- grep(dmiHAD$Season[1], oni$Season)
    lr <- grep(dmiHAD$Season[nrow(dmiHAD)], oni$Season)
    dmiHAD <- cbind(oni[fr:lr,1:3],dmiHAD[,2:13])
    colnames(dmiHAD)[4:15] <- colnames(oni)[4:15]
    
    # Combine aoi data sets
    return(list(ONI = oni, MEI = mei, DMI = dmi, DMIHAD = dmiHAD))
  }
  
}