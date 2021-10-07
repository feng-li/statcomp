library(shiny)

# Define UI 
fluidPage(
  
  # Application title
  titlePanel("Please choose a mean and sd"),
  
  # Sidebar with controls 
  sidebarLayout(
    sidebarPanel(
      numericInput("mean", 
                   label = h3("Mean"), 
                   value = 1),
      
      selectInput("sd", label = h3("Standard Deviation"), 
                  choices = list("1" = 1, "2" = 2,
                                 "3" = 3), selected = 1)
  
    ),
    

    mainPanel(
      h2("This is the density plot of normal distribution."),
      h3("Note how it changes with mean and sd."),
      plotOutput("normPlot")
    )
  )
)
