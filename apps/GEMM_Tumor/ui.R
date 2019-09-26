library(shiny)

shinyUI(fluidPage(
  titlePanel("Gene expression from RNA-seq of GEMM ESR1 mutant tumors"),
  
  sidebarLayout(
    sidebarPanel(
      
      textInput("gene", label= "Enter Gene Symbol in CAPITAL letters", value = "ESR1"),
      strong("RNA-seq data from 9 ESR1 WT tumors and 8 ESR1 Y541S tumors"),
      helpText("RNAseq was done on ESR1 WT/Mutant tumors samples uisng TrueSeq Library protocol 
               llumina TruSeq Stranded Total RNA Kit , NovaSeq SP300 flow-cell with 101 paired-dnd reds.
               Among all the tumors, WT 9 and 10; Mut 7 and 8; Mut 9 and 10 were from same mice"),
      tableOutput("Table"),

    
    width=4),
    
    mainPanel(
      fluidRow(
        column(5,plotOutput("ER"))),
      fluidRow(
        column(5,plotOutput("ER2")))
      
    ))
    
  
)
)