load("LTED2WT_log2fc_adjp.rda")
load("LTEDvs1week_log2fc_adjp.rda")
load("res_SUM44_MM134_all.rda")
load("SUM44_MM134_LTED_RNAseq_gene_TPM_untrimmed.rda")
shinyServer(function(input, output) {
  queryExpression <- reactive({
    queryExpression_LTEDvsWT <- LTED2WT_log2fc_p[LTED2WT_log2fc_p$gene_symbol %in% input$geneExpressionQuery,]
    queryExpression_LTEDvs1week <- LTEDvs1week_log2fc_p[LTEDvs1week_log2fc_p$gene_symbol %in% input$geneExpressionQuery,]
    queryExpression<-rbind(queryExpression_LTEDvsWT,queryExpression_LTEDvs1week)
    queryExpression$Type<-c("LTEDvsWT","LTEDvs1week")
    queryExpression<-queryExpression[,c(12,9,10,1,2,3,4,5,6,7,8)]
    return(queryExpression)
  })
 output$expression<-renderDataTable({
    return(queryExpression())
  }, options = list(lengthMenu = c(10,20), pageLength=40)
  )
output$LTEDgeneFoldchange<-renderDataTable({
   gene<-input$geneExpressionQuery
   df<-as.data.frame(matrix(NA,2,6))
   colnames(df)<-c("MM134L/A","MM134L/B","MM134L/D","MM134L/E","SUM44L/A","SUM44L/B")
   rownames(df)<-c("Log2FoldChange","adjusted p-value")
   df[1,]<-res_SUM44_MM134_all[which(res_SUM44_MM134_all$gene_symbol==gene),c(1,3,5,7,9,11)]
   df[2,]<-res_SUM44_MM134_all[which(res_SUM44_MM134_all$gene_symbol==gene),c(2,4,6,8,10,12)]
   df<-t(df)
   df<-cbind(c("MM134L/A","MM134L/B","MM134L/D","MM134L/E","SUM44L/A","SUM44L/B"),df)
   return(df)
 })
output$expressionPlot<-renderPlot({
  gene<-input$geneExpressionQuery
  df<-t(res_SUM44_MM134_all[which(res_SUM44_MM134_all$gene_symbol==gene),c(1,3,5,7,9,11)])
  df$expression<-2^df[,1]
  barplot(df$expression, main =paste(gene,"Foldchange LTED vs parental",sep = " ") , col = c("red","red","red","red","blue","blue"),names =c("134L/A","134L/B","134L/D","134L/E","44L/A","44L/B"),cex.names = 1.2,las=2,ylim = c(0,1.2*max(df$expression)))
  abline(h=1, col="blue", lwd=2)
})  
output$TPMplot <- renderPlot({
  gene<-input$geneExpressionQuery
  df<-t(LTED.cell.line.TPM_gene_untrimmed[which(LTED.cell.line.TPM_gene_untrimmed$gene_symbol==gene),c(-16,-17,-18,-28)])
  barplot(df[,1], main = "TPM of parental and LTED cells", col = c("black","black","black","red","red","red","red","red","red","red","red","red","red","red","red","grey","grey","grey","blue","blue","blue","blue","blue","blue"),names =c("134P.1","134P.2","134P.3","134LA.1","134LA.2","134LA.3","134LB.1","134LB.2","134LB.3","134LD.1","134LD.2","134LD.3","134LE.1","134LE.2","134LE.3","44P.1","44P.2","44P.3","44LA.1","44LA.2","44LA.3","44LB.1","44LB.2","44LB.3"),cex.names = 1.2,las=2,ylim=c(0,1.1*max(df[,1])))
  abline(h=1, col="blue", lwd=2) 
})
output$LTEDvsWTPlot <- renderPlot({
  gene<-input$geneExpressionQuery
  Genedata<-as.data.frame(t(LTED2WT_log2fc_p[which(LTED2WT_log2fc_p$gene_symbol %in% input$geneExpressionQuery),c(9,1,3,5,7)]))
  colnames(Genedata)<-"expression"
  Genedata$expression<-2^Genedata$expression
  barplot(Genedata$expression, main =paste(gene,"expression change LTED vs WT",sep=" "),ylab="Foldchange",col = c("red","blue","blue","grey","grey"),names =c("SUM44","HCC1428","MCF7","T47D","ZR751"),cex.lab=1.1,cex.names = 1.1,las=2,ylim=c(0,1.2*max(Genedata$expression)))
  abline(h=1, col="black", lwd=2)
})
output$LTEDvs1weekPlot <- renderPlot({
  gene<-input$geneExpressionQuery
  Genedata<-as.data.frame(t(LTEDvs1week_log2fc_p[which(LTEDvs1week_log2fc_p$gene_symbol %in% input$geneExpressionQuery),c(9,1,3,5,7)]))
  colnames(Genedata)<-"expression"
  Genedata$expression<-2^Genedata$expression
  barplot(Genedata$expression, main =paste(gene,"expression change LTED vs 1 week HD",sep=" "),ylab="Foldchange",col = c("red","blue","blue","grey","grey"),names =c("SUM44","HCC1428","MCF7","T47D","ZR751"),cex.lab=1.1,cex.names = 1.1,las=2,ylim=c(0,1.2*max(Genedata$expression)))
  abline(h=1, col="black", lwd=2)
})
})
