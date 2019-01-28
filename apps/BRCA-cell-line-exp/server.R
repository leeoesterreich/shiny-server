load("Neel_cell_line_per_gene_normalized_log2_FPKM_Allgene.rda")
load("Neel_cell_line_subtype.rda")
load("Cell_line_zGARP_Allgene.rda")
shinyServer(function(input, output) {
  querydata <- reactive({
    queryExpression <- RNAseq_normalized[rownames(RNAseq_normalized) %in% input$geneExpressionQuery,]
    queryExpression <- queryExpression[,-c(1,2,3)]
    queryzGARP<-zGARP[rownames(zGARP) %in% input$geneExpressionQuery,]
    queryzGARP<-queryzGARP[,-c(1,2,3)]
    queryzGARPdf<-t(queryzGARP)
    queryExpressiondf<-t(queryExpression)
    querydata<-as.data.frame(matrix(NA,81,2))
    rownames(querydata)<-rownames(queryExpressiondf)
    colnames(querydata)<-c("Expression","zGARP")
    querydata[,1]<-queryExpressiondf[rownames(querydata),1]
    querydata[rownames(queryzGARPdf),2]<-queryzGARPdf[,1]
    return(querydata)
  })

  final.data <- reactive({
    querydata<- querydata()
    dataordered<-querydata[order(querydata[,1],na.last = F),]
    final.data<-cbind(rownames(dataordered),dataordered)
    colnames(final.data)<-c("Cell_line","Log2FPKM","zGARP")
    return(final.data)
  })
  output$expression <- renderDataTable({
    return(final.data())
  }, options = list(lengthMenu = c(50,100), pageLength=100)
  )
  
  observe({
    if(input$plotCommand == TRUE) {
      output$expressionPlot<-renderPlot({
        gene<-input$geneExpressionQuery
        Genedata<-as.data.frame(t(data.matrix(RNAseq_normalized[gene,-c(1,2,3)])))
        Genedata$histology<-cell_line_subtype[rownames(Genedata),"histology"]
        df<-as.data.frame(Genedata)
        colnames(df)<-c("gene", "histology")
        p<-wilcox.test(df[which(df$histology=="ER+ILC"),"gene"], df[which(df$histology=="Others"),"gene"], paired = F, alternative = "two.sided")$p.value
        p<-paste("Mann-Whitney U p-value=", signif(p, digits=2))
        boxplot(gene ~ histology, data = df, names=c("ER+ILC","Others"), outpch=NA, xlab=p)
        stripchart(gene~ histology, data = df, vertical = TRUE, method = "jitter",
                   pch = 21, col = c("red", "blue"),add = TRUE) 
        title(paste(gene, ": All cell lines ER+ILC vs. Others"))
      })
    }
  })
  
  
  observe({
    if(input$barPlot == TRUE) {
      output$expressionPlot <- renderPlot({
        gene<-input$geneExpressionQuery
        Genedata<-as.data.frame(t(data.matrix(RNAseq_normalized[gene,-c(1,2,3)])))
        Genedata$expression<-2^Genedata[,1]
        Genedata_cell_line<-Genedata[c("MDAMB134","SUM44","MDAMB330","T47D","MCF7","BT474","MDAMB231","MCF10A","MCF12A"),]
        barplot(Genedata_cell_line$expression, main = "Normalized gene expression (FPKM)", col = "red",names =c("MM134","SUM44","MM330","T47D","MCF7","BT474","MM231","MCF10A","MCF12A"),cex.names = 1.2,las=2)
      })
    }
  })
  
  observe({
    if(input$plotzGARP == TRUE) {
      output$expressionPlot <- renderPlot({
        gene<-input$geneExpressionQuery
        Genedata<-as.data.frame(t(data.matrix(zGARP[gene,-c(1,2,3)])))
        Genedata_cell_line<-Genedata[c("MDAMB134","SUM44","MDAMB330","T47D","MCF7","BT474","MDAMB231","MCF10A","MCF12A"),]
        barplot(Genedata_cell_line, main = "zGARP Score", col = "red",names =c("MM134","SUM44","MM330","T47D","MCF7","BT474","MM231","MCF10A","MCF12A"),cex.names = 1.2,las=2)
      })
    }
  })
})



