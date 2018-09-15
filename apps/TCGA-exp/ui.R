#library(shiny)
#load("primary_tpm.Rda")
#genelist<-as.list(colnames(df[,1:23368]))  # 1:23368

shinyUI(fluidPage(
  titlePanel("  TCGA Gene Expression (TPM) Comparison"),
  
  sidebarPanel(
    #selectInput(inputId="gene", label = "Gene Query", choices = c(genelist), selected = "ESR1"),
    textInput("gene", "Gene Query", value = "ESR1"),
    selectInput(inputId="dataset", 
                label = "Study Population",
                choices = list("Breast Cancer", 
                               "ER+",
                               "ER-",
                               "ILC",
                               "IDC",
                               "LumA",
                               "LumB",
                               "Her2",
                               "Basal"),
                selected = "Breast Cancer"),
    selectInput(inputId="group",
                label = "Comparison Group",
                choices = list("ER+ vs ER-","ILC vs IDC","LumA vs LumB"),
                selected = "ER+ vs ER-"),
    
    verbatimTextOutput("state"),
    verbatimTextOutput("stats1"),
    verbatimTextOutput("stats2"),
    verbatimTextOutput("stats3"),
    verbatimTextOutput("stats4")
  ),
  mainPanel(
    fluidRow(
      splitLayout(cellWidths = c("65%", "35%"), plotOutput("plot1"), plotOutput("plot2"))
    ),
    fluidRow(
      splitLayout(cellWidths = c("15%", "70%", "15%"),plotOutput("a"),plotOutput("plot3"),plotOutput("b") )
    ),
    
    dataTableOutput("view1")
  )
)
)