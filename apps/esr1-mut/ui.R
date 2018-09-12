library(shiny)

shinyUI(fluidPage(
  titlePanel("Plot TPM"),
  
  sidebarLayout(
    sidebarPanel(
      
      textInput("gene", label= "Enter Gene Symbol in all CAPITAL letters", value = "ESR1")
      
      
    ),
    
    mainPanel(plotOutput("plotTPM_MCF7"), plotOutput("plotTPM_T47D"))
    #mainPanel(
    #         fluidRow(column(width=6, plotOutput("plotTPM_MCF7")),
     #          column(width=6,plotOutput("plotTPM_T47D"))
    #         )
    #  )
    )
  )
  )