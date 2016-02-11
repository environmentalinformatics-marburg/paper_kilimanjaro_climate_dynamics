library(shiny)

load("trend_analysis.RData")
values <- colnames(rf[colnames(rf) != "YEAR"])

str(rf)
for(i in c(9, 15:34)){
rf[, i] <- as.factor(rf[, i] )
}


shinyServer(function(input, output) {
  
  output$dep_prm <- renderUI({
    selectInput("dep_prm", "Dependent parameters:", 
                values, selected = "P_MEAN_NEW")
  })
  
  output$ind_prm <- renderUI({
    checkboxGroupInput("ind_prm", "Independent parameters:", values)
  })
  
  frml <- reactive({
    if(length(input$ind_prm) == 0){
      as.formula(paste(input$dep_prm," ~ time(", input$dep_prm, ")"))
    } else if(length(input$ind_prm) == 1){
      as.formula(paste0(input$dep_prm," ~ time(", input$dep_prm, ")", 
                        paste0("+ ", input$ind_prm, " * time(", input$dep_prm, ")")))
    } else {
      as.formula(paste0(input$dep_prm," ~ time(", input$dep_prm, ")", 
                        paste0("+ ", input$ind_prm, collapse = paste0(" * time(", input$dep_prm, ")"))))
    }
  })
  
  rflm <- reactive({lm(frml(), data = rf)})
                       
  
  output$frml <- renderPrint(frml())
  output$lm <- renderPrint({
    summary(rflm())
  })
  output$anova <- renderPrint({
    summary(rflm())
  })
})



