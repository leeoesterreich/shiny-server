library(shiny)
library(dqshiny)
load("autoFillNames.Rda")

shinyUI(fluidPage(

  titlePanel("BrCa Cell Line Expression"),

  sidebarLayout(
    
    sidebarPanel(
      autocomplete_input("geneEntry", "Gene Query", value = "ERBB2_ENSG00000141736", autoFillNames, placeholder = "Type Gene Name"),
      p("Data generated from RNA-seq fastqs using salmon (PMID: 28263959, v0.12.0, 31-kmer quasi-mapping, gcBias, seqBias, validateMappings). Fastqs from either (1) CCLE (PMID: 22460905), (2) Marcotte et al Cell 2016 (PMID: 26771497), (3) courtesy of Dr. Joe Gray, OHSU. Mapping reference: Homo_sapiens.GRCh38.82. Expression represented as log2-transformed TMM-normalized CPMs. PAM50s called with genefu (PMID: 26607490)"),
      checkboxInput("expPlot", "Plot cell line data", FALSE),
      checkboxInput("distributionPlot", "Plot total distribution", FALSE),
      checkboxInput("pam50plot", "Plot PAM50 expression distribution", FALSE),
      plotOutput("expressionPlot"), width = 4
        ),

   mainPanel(
      dataTableOutput("plotData")
    )
   )
))