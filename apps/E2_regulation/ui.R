library(shiny)

shinyUI(fluidPage(
  titlePanel("Gene Expression in ER+ Breast Cancer Cell Lines -/+ Estrogen"),
  
  sidebarLayout(
    sidebarPanel(
      
      textInput("gene", label= "Enter Gene Symbol in CAPITAL letters", value = "ESR1"),

      helpText("All data were extracted from the following resources:"),
      strong("Disclaimer to users: Inter-data-sets normalization was not performed. Please do not compare magnitude between different data sets."),        
      tableOutput("Table"),
    width=3),
    
    mainPanel(
      fluidRow(
        column(3,plotOutput("MCF7_1")), 
        column(3,plotOutput("MCF7_2")), 
        column(3,plotOutput("MCF7_3")),
        column(3,plotOutput("MCF7_4"))),
      fluidRow(
        column(3,plotOutput("T47D_1")),
        column(3,plotOutput("T47D_2")), 
        column(3,plotOutput("T47D_3")),
        column(3,plotOutput("T47D_4"))),
      fluidRow(
        column(3,plotOutput("ZR")),
        column(3,plotOutput("BT474")), 
        column(3,plotOutput("MM134")),
        column(3,plotOutput("SUM44")))
    
        
      ))
    
  
)
)
