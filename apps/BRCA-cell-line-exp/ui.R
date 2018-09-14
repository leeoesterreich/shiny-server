library(shiny)

shinyUI(fluidPage(

  titlePanel("BrCa Cell Line Normalized Log2 Expression"),

  sidebarLayout(
    
    sidebarPanel(
      textInput("geneExpressionQuery", "Gene Query", value = "ESR1"),
      p("Data from Marcotte et al 2016. Per gene log2 FPKM was Quantile-Normalized, and a floor value set at log2(0.1) FPKM.Per gene zGARP score was served as gene essentiality scores. 3 ER+ILC cells lines MDAMB134, MDAMB330, SUM44 are in the dataset"),
      checkboxInput("plotCommand", "Plot Normalized log2 FPKM in ER+ILC cell lines and other cell lines?", FALSE),
      checkboxInput("barPlot", "Plot normalized FPKM data of SUM44, MM134, MM330, T47D, MCF7, BT474, MDAMB231, MCF10A, MCF12A?", FALSE),
      checkboxInput("plotzGARP", "Plot gene zGARP score of SUM44, MM134, MM330, T47D, MCF7, BT474, MDAMB231, MCF10A, MCF12A?", FALSE),  
      plotOutput("expressionPlot")
      ),
   mainPanel(
      dataTableOutput("expression")
    )
   )
))