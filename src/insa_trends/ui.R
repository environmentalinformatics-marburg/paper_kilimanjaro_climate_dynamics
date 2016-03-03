library(shiny)



# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Insa's trends"),
  
  # Sidebar and mainbar
  sidebarLayout(
    position = "right",
    sidebarPanel(
      uiOutput("dataset"),
      uiOutput("dep_prm"),
      uiOutput("ind_prm")
    ),
    
    # Show plots
    mainPanel(
      textOutput("frml"),
      verbatimTextOutput ("lm"),
      verbatimTextOutput("sw"),
      plotOutput("plot"),
      plotOutput("qqplot"),
      verbatimTextOutput ("anova")
      # plotOutput("distPlot")
    )
  )
))

# # Rainfall annual trends
# for (x in unique(rfa$Parameter)){
#   print(x)
#   print(summary(lm(Rainfall ~ time(Year), data = rfa[rfa$Parameter == x,])))
# }
# print(summary(rq(Rainfall ~ time(Year), tau = seq(0.8, 1, 0.05),
#                  data = rfa[rfa$Parameter == "P_MEAN_NEW",])))
# 
# # Time series stationarity
# kpss.test(rf$P_MEAN_NEW)
# kpss.test(rf$P_MEAN_NEW, null = "Trend")
# 
# 
# # Rainfall frequency
# k = kernel("daniell", c(3,3))
# rfspec <- spec.pgram(rf$P_MEAN_NEW, kernel = k, taper = 0, plot = FALSE, detrend = TRUE)
# 
# rfspecdf <- data.frame(freq = rfspec$freq, spec = rfspec$spec)
# ticks <- c(30, 20, 10, 5, 3, 1, 0.5, 0.25)
# ticks <- c(30, 20, seq(12, 1, -1), seq(1, 0, -0.25))
# labels <- as.character(ticks)
# 
# breaks <- 1/(ticks * 12)
# 
# ggplot(rfspecdf, aes(x = freq, y = spec)) + 
#   geom_line() + 
#   scale_x_log10("Period (years)", breaks = breaks, labels = labels) + 
#   scale_y_log10()
# 
# 
# rfwave <- morlet(y1 = rf$P_MEAN_NEW, x1 = time(rf$YEAR), p2 = 8, dj = 0.1, siglvl = 0.95)
# rfwave$period <- rfwave$period/12
# levels <- quantile(rfwave$Power, c(0, 0.25, 0.5, 0.75, 0.95, 1))
# wavelet.plot(rfwave, wavelet.levels = levels, crn.ylim = c(22.5, 30))
# 
# rfwave_avg <- data.frame(power = apply(wave.out$Power, 2, mean), period = (wave.out$period))
# ggplot(rfwave_avg, aes(x = period, y = power)) + 
#   geom_line() + 
#   scale_x_continuous(breaks = seq(25)) + 
#   scale_y_continuous()
# 
# # Rainfall frequency de-seasoned
# k = kernel("daniell", c(3,3))
# rfspec <- spec.pgram(rf$P_MEAN_NEW_ds, kernel = k, taper = 0, plot = FALSE, detrend = TRUE)
# 
# rfspecdf <- data.frame(freq = rfspec$freq, spec = rfspec$spec)
# ticks <- c(30, 20, 10, 5, 3, 1, 0.5, 0.25)
# ticks <- c(30, 20, seq(12, 1, -1), seq(1, 0, -0.25))
# labels <- as.character(ticks)
# 
# breaks <- 1/(ticks * 12)
# 
# ggplot(rfspecdf, aes(x = freq, y = spec)) + 
#   geom_line() + 
#   scale_x_log10("Period (years)", breaks = breaks, labels = labels) + 
#   scale_y_log10()
# 
# 
# rfwave <- morlet(y1 = rf$P_MEAN_NEW_ds, x1 = time(rf$YEAR), p2 = 8, dj = 0.1, siglvl = 0.95)
# rfwave$period <- rfwave$period/12
# levels <- quantile(rfwave$Power, c(0, 0.25, 0.5, 0.75, 0.95, 1))
# wavelet.plot(rfwave, wavelet.levels = levels, crn.ylim = c(22.5, 30))
# 
# rfwave_avg <- data.frame(power = apply(wave.out$Power, 2, mean), period = (wave.out$period))
# ggplot(rfwave_avg, aes(x = period, y = power)) + 
#   geom_line() + 
#   scale_x_continuous(breaks = seq(25)) + 
#   scale_y_continuous()
# 
# 
# 
# # Yearly rainfall sums
# ggplot(rf, aes(x = mnth, y = P_MEAN_NEW, colour = as.factor(Season))) +
#   geom_bar(stat = "identity", position = "dodge") +
#   xlab("Year") +
#   ylab("Precipitation")
# 
# # Rainfall trends
# rflm <- lm(P_MEAN_NEW_ds ~  time(YEAR), data = rf[rf$mnth == 4,])
# summary(rflm)
# anova(rflm)
# 
# plot(rflm)
# 
# rflm <- lm(P_MEAN_NEW_SQRT2 ~  time(YEAR) + short_rains * time(YEAR), data = rf)
# summary(rflm)
# anova(rflm)
# 
# rflm <- lm(P_MEAN_NEW_SQRT2 ~  time(YEAR) + long_rains * time(YEAR), data = rf)
# summary(rflm)
# anova(rflm)
# 
# rflm <- lm(P_MEAN_NEW_SQRT2 ~  time(YEAR) + long_rains * time(YEAR) + short_rains * time(YEAR), data = rf)
# summary(rflm)
# anova(rflm)
# 
# 
# 
# rflm <- lm(P_MEAN_NEW_ds ~  time(YEAR) + as.factor(mnth) * time(YEAR), data = rf)
# summary(rflm)
# anova(rflm)
# 
# rflm <- lm(P_MEAN_NEW_LOG ~  time(YEAR) + drywet * time(YEAR), data = rf)
# summary(rflm)
# anova(rflm)
# 
# 
# 
# 
# 
# 
# 
# 
# rflm <- lm(P_MEAN_NEW_ds ~  time(YEAR) + long_rains * time(YEAR) + short_rains * time(YEAR) +
#              M1 * time(YEAR) +
#              M3 * time(YEAR) + M4 * time(YEAR) +  M5 * time(YEAR) + 
#              M6 * time(YEAR) + M7 * time(YEAR) + M8 * time(YEAR) + 
#              M9 * time(YEAR) + M10 * time(YEAR) + M11 * time(YEAR) +
#              M12 * time(YEAR), data = rf)
# summary(rflm)
# anova(rflm)
# 
# rflm <- lm(P_MEAN_NEW_ds ~  time(YEAR) + long_rains * time(YEAR) + short_rains * time(YEAR) +
#              M1 * time(YEAR) + M2 * time(YEAR) +
#              M3 * time(YEAR) + M4 * time(YEAR) +  M5 * time(YEAR) + 
#              M6 * time(YEAR) + M7 * time(YEAR) + M8 * time(YEAR) + 
#              M9 * time(YEAR) + M10 * time(YEAR) + M11 * time(YEAR) +
#              M12 * time(YEAR) + 
#              IOD * time(YEAR) + ENSO * time(YEAR), data = rf)
# summary(rflm)
# anova(rflm)
# 
# rflm <- lm(P_MEAN_NEW_LOG ~  time(YEAR) + long_rains * time(YEAR) + short_rains * time(YEAR) +
#              M1 * time(YEAR) + M2 * time(YEAR) +
#              M3 * time(YEAR) + M4 * time(YEAR) +  M5 * time(YEAR) + 
#              M6 * time(YEAR) + M7 * time(YEAR) + M8 * time(YEAR) + 
#              M9 * time(YEAR) + M10 * time(YEAR) + M11 * time(YEAR) +
#              M12 * time(YEAR) + 
#              IOD * time(YEAR) + ENSO * time(YEAR), data = rf)
# summary(rflm)
# anova(rflm)
# 
# rflm <- lm(P_MEAN_NEW ~  time(YEAR) + long_rains * time(YEAR) + short_rains * time(YEAR) +
#              M1 * time(YEAR) + M2 * time(YEAR) +
#              M3 * time(YEAR) + M4 * time(YEAR) +  M5 * time(YEAR) + 
#              M6 * time(YEAR) + M7 * time(YEAR) + M8 * time(YEAR) + 
#              M9 * time(YEAR) + M10 * time(YEAR) + M11 * time(YEAR) +
#              M12 * time(YEAR) + 
#              IOD * time(YEAR) + ElNino * time(LaNina), data = rf)
# summary(rflm)
# anova(rflm)
