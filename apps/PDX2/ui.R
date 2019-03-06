library(shiny)

shinyUI(fluidPage(
  titlePanel("Gene expression from RNA-seq of Patient Derived Xenograft models (HCI001-HCI019)"),
  
  sidebarLayout(
    sidebarPanel(
      
      textInput("gene", label= "Enter Gene Symbol in CAPITAL letters", value = "ESR1"),
      strong("RNA-seq data from HCI001-HCI019 PDX models were extracted from GSE 113476:"),
      helpText("RNAseq was done on Breast cancer PDX samples uisng Library protocol 
               =llumina TruSeq Stranded Total RNA Kit with Ribo-Zero Gold , HiSeq 125 Cycle Paired-End Sequencing v4"),
      tableOutput("Table"),

    
    width=3),
    
    mainPanel(
      fluidRow(
        column(7,plotOutput("Types"))),
      fluidRow(
        column(7,plotOutput("ER")))
        
      ))
    
  
)
)