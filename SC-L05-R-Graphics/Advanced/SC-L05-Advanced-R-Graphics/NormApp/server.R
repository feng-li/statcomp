library(shiny)


# Define server logic required to plot density of normal distribution
function(input, output) {
  

  output$normPlot <- renderPlot({
    plot(density(rnorm(1000, mean = as.numeric(input$mean), 
               sd = as.numeric(input$sd))), main = 'normal density',
         xlim = c(-10, 10))
  })
}