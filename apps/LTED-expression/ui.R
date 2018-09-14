library(shiny)

shinyUI(fluidPage(

  titlePanel("LTED cell lines expression foldchange"),

  sidebarLayout(
    
    sidebarPanel(
      textInput("geneExpressionQuery", "Gene Query", value = "ESR1"),
      p("Data from Simigdala et al, BCR 2016. WT parental cells are cutured in RPMI pheno red free +10%FBS+1nM E2. Parental cells undergone 1 week hormone deprivation, WT parental, and LTED cells were performed microarray.For lab LTED, parental cell lines were hormone deprived for 3 days before RNASeq. Salmon and DESeq2 were used to map reads to genes and select the differentially expressed genes. SUM44F is used as SUM44 parental cells"),
      plotOutput("LTEDvsWTPlot"),
      plotOutput("LTEDvs1weekPlot")
      ),
   mainPanel(
      dataTableOutput("expression"),dataTableOutput("LTEDgeneFoldchange"),plotOutput("expressionPlot"),plotOutput("TPMplot")
    )
   )
))