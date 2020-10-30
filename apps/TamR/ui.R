library(shiny)

shinyUI(fluidPage(
  titlePanel("Gene Expression in Tamoxifen or Fulvestrant Resistant Cell Models"),
  
  sidebarLayout(
    sidebarPanel(
      
      textInput("gene", label= "Enter Gene Symbol in CAPITAL letters", value = "ESR1"),

      helpText("All data were extracted from the following resources:"),
      strong("Disclaimer to users: Inter-data-sets normalization was not performmed. Please do not compare magnitude between data sets."),        
      tableOutput("Table"),
    width=3),
    
    mainPanel(
      fluidRow(
        column(3,plotOutput("DF1")), 
        column(3,plotOutput("DF2")), 
        column(3,plotOutput("DF4")),
        column(3,plotOutput("DF5"))),
      fluidRow(
        column(3,plotOutput("DF6")),
        column(3,plotOutput("DF7_1")), 
        column(3,plotOutput("DF7_2")),
        column(3,plotOutput("DF7_3"))),
      fluidRow(
        column(3,plotOutput("DF7_4")),
        column(3,plotOutput("DF8")), 
        column(3,plotOutput("DF11")),
        column(3,plotOutput("DF12")))
    
        
      ))
    
  
)
)