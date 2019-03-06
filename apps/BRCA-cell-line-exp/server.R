load("brcaLines.tmmNorm.log2cpm.shiny.Rda")
load("autoFillNames.Rda")
library(ggplot2)
library(ggbeeswarm)

shinyServer(function(input, output) {
  
  plotData <- reactive({  
  plotData <- brcaLines.tmmNorm.log2cpm.shiny[,c(input$geneEntry,"PAM50")]
  plotData$cellLine <- rownames(plotData)
  plotData.ordered <- plotData[order(-plotData[,1]),]
  plotData.ordered$cellLine <- factor(plotData.ordered$cellLine, levels = plotData.ordered$cellLine)
  plotData <- plotData.ordered
    return(plotData)
    })

  output$plotData <- renderDataTable({
    return(plotData())
    }, options = list(lengthMenu = c(50,100), pageLength=100)
    )

  observe({
    if(input$expPlot == TRUE) {
    output$expressionPlot <- renderPlot({

  plotData <- brcaLines.tmmNorm.log2cpm.shiny[,c(input$geneEntry,"PAM50")]
  plotData$cellLine <- rownames(plotData)
  plotData.ordered <- plotData[order(-plotData[,1]),]
  plotData.ordered$cellLine <- factor(plotData.ordered$cellLine, levels = plotData.ordered$cellLine)
  plotData <- plotData.ordered
 
  ggplot(plotData, aes_string(x = "cellLine", y = input$geneEntry)) +
  geom_point(size=4, alpha = 0.7, aes_string(col = "PAM50")) +
  ggtitle(paste("BrCa Lines", input$geneEntry, "Expression")) + 
  theme(
    panel.border = element_rect(fill = NA,color = "black", size = 2),
    panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size=10, face = "bold"), 
    axis.text.x = element_text(size = 6, angle = 90, vjust=1),
    axis.text.y = element_text(size = 8),
    panel.background = element_rect(fill = "gray92")) + 
    ylab(paste(input$geneEntry, 'Expression (normlog2CPM)', sep = " "))
      })
    }
  })


observe({

    if(input$distributionPlot == TRUE) {
    output$expressionPlot <- renderPlot({

  plotData <- brcaLines.tmmNorm.log2cpm.shiny[,c(input$geneEntry,"PAM50")]

   data <- as.numeric(plotData[,1])
   boxplot(data, outpch=NA, frame = T, main = "Expression Distribution", border = "blue3")
   return(stripchart(data, vertical = TRUE, method = "jitter", pch = 21, col = "red", add = TRUE, log="y"))


      })
    }
  })

  observe({


    if(input$pam50plot == TRUE) {
    output$expressionPlot <- renderPlot({
      
        plotData <- brcaLines.tmmNorm.log2cpm.shiny[,c(input$geneEntry,"PAM50")]
        plotData$cellLine <- rownames(plotData)
        plotData.ordered <- plotData[order(-plotData[,1]),]
        plotData.ordered$cellLine <- factor(plotData.ordered$cellLine, levels = plotData.ordered$cellLine)
        plotData <- plotData.ordered

      ggplot(plotData, aes_string(y = input$geneEntry, x = "PAM50", fill = "PAM50" )) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
      geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.3, alpha = 0.5) + 
      ggtitle(paste("BrCa Lines", input$geneEntry, "PAM50 Expression")) +
        theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 0, face = "bold"), 
          axis.text.x = element_text(size = 8, face = "bold", angle = 90),
          axis.title.y = element_text(size=8, face = "bold"),
          panel.background = element_rect(fill = "gray92"),
          legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
       xlab('PAM50') +
       ylab(paste(input$geneEntry, 'Expression (normlog2CPM)', sep = " "))

      })
    }
  })


})
