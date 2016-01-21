#### Set working environment ###################################################

library(ggplot2)

setwd("D:/active/paper_kilimanjaro_climate_dynamics")
rf <- read.table("metoffice_1973-2013.csv", header = TRUE, sep = ";", dec = ".")

rf$YEAR <- as.Date(as.character(rf$YEAR))

head(rf)

rfa <- lapply(seq(2, length(colnames(rf))), function(x){
  ret <- aggregate(rf[, x], by = list(substr(rf$YEAR, 1, 4)), FUN = "sum")
  ret$Group.1 <- as.Date(ret$Group.1, format = "%Y")
  colnames(ret) <- c("Year", "Rainfall")
  ret$Station <- colnames(rf)[x]
  return(ret)
})
rfa <- do.call("rbind", rfa)


for (x in unique(rfa$Station)){
  print(x)
  print(summary(lm(Rainfall ~ time(Year), data = rfa[rfa$Station == x,])))
}

