library(shiny)
# library(dplR)
library(ggplot2)
# library(multitaper)
# library(quantreg)
# library(tseries)

if(Sys.info()["user"] == "shiny"){
  in_path <- "/srv/shiny-server/moc/temp/insa_trends/"
  load("trend_analysis.RData")
} else {
  in_path <- "D:/active/paper_kilimanjaro_climate_dynamics/src/insa_trends/"
  load(paste0(in_path, "trend_analysis.RData"))
}




shinyServer(function(input, output, session) {
  
  output$dataset <- renderUI({
    selectInput("dataset", "Dataset:", 
                c("rf", "rfa"), selected = "rf")
  })
  
  output$dep_prm <- renderUI({
    selectInput("dep_prm", "Dependent parameters:", 
                values(), selected = "P_MEAN_NEW")
  })
  
  output$ind_prm <- renderUI({
    checkboxGroupInput("ind_prm", "Independent parameters:", values())
  })
  
  dfin <- reactive({
    df <- switch(input$dataset,
           "rf" = rf,
           "rfa" = rfa)
    if(input$dataset == "rf"){
      for(i in c(9, 17:36)){
        df[, i] <- as.factor(df[, i] )
      }
    } else {
      df <- df[df$Parameter == "P_MEAN_NEW", ]
    }
    return(df)
  })

  values <- reactive({
    df <- dfin()
    colnames(df[colnames(df) != "Year"])})
  
  frml <- reactive({
    if(length(input$ind_prm) == 0){
      as.formula(paste(input$dep_prm," ~ time(", input$dep_prm, ")"))
    } else if(length(input$ind_prm) == 1){
      as.formula(paste0(input$dep_prm," ~ time(", input$dep_prm, ")", 
                        paste0("+ ", input$ind_prm, " * time(", input$dep_prm, ")")))
    } else {
      as.formula(paste0(input$dep_prm," ~ time(", input$dep_prm, ") +", 
                        paste0(input$ind_prm, collapse = paste0(" * time(", input$dep_prm, ") + ")),
                        paste0(" * time(", input$dep_prm, ")")))
    }
  })
  
  
  dflm <- reactive({
    df <- dfin()    
    lm(frml(), data = df)
  })
  
  lmplot <- reactive({
    df <- dfin()
    df$pred <- predict(dflm())
    df$time <- time(df$Year)
    if(length(input$ind_prm) == 0){
      ggplot(df, aes_string(x = "time", y = input$dep_prm)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_line(aes(y = pred)) + 
        xlab("Year") +
        ylab("Precipitation")
    } else {
      ggplot(df, aes_string(x = "time", y = input$dep_prm), environment = environment() ) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_line(aes(y = pred, color = interaction(df[, input$ind_prm]))) + 
        xlab("Year") +
        ylab("Precipitation")
    }
  })
  
  
  qqplot <- reactive({
    vec <- rstandard(dflm())
    y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
    x <- qnorm(c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
    d <- data.frame(resids = vec)
    ggplot(d, aes(sample = resids)) + stat_qq() + geom_abline(slope = slope, intercept = int)
  })
  
  
  output$frml <- renderPrint(frml())
  
  output$lm <- renderPrint({
    summary(dflm())
  })
  
  output$sw <- renderPrint({
    shapiro.test(residuals(dflm()))
  })
  
  output$plot <- renderPlot({
    lmplot()
  })
  
  output$qqplot <- renderPlot({
    qqplot()
  })
  
  output$anova <- renderPrint({
    summary(dflm())
  })
})



