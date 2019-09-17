library(shiny)
library(ggplot2)

GEMM_CPM<-read.csv("GEMM_CPM.csv")
group<-read.csv("Group.csv",header=T)


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  output$Table<- renderTable(group)
   
 
 output$ER<- renderPlot({
   gene1<-input$gene
   rownames(GEMM_CPM)=make.unique(as.character(GEMM_CPM$Human))
   GEMM_CPM=GEMM_CPM[,2:18]
   GEMM_CPM1<-t(GEMM_CPM)
   class(GEMM_CPM1)
   GEMM_CPM1=as.data.frame(GEMM_CPM1)
   GEMM_CPM1$ER<-group$ER
   df<-GEMM_CPM1[,c("ER",gene1)]
   colnames(df)<-c("ER","gene1")
   
   p2<-wilcox.test(df[which(df$ER=="WT"),"gene1"], df[which(df$ER=="Mutant"),"gene1"], paired = F, alternative = "two.sided")$p.value
   p2<-paste("Mann-Whitney U p-value (WT vs Mut)=", signif(p2, digits=2))
   
   boxplot(gene1 ~ ER, data=df,outpch=NA,ylab="TMMLog2(CPM+1)", par(cex.lab=1.5),par(cex.axis=1.5),par(cex.main=2),xlab=p2)
   stripchart(gene1 ~ ER, data = df, vertical = TRUE, method = "jitter",
              pch = 20, col = c("red","black"),add = TRUE) 
   title(paste(gene1, "Expression~ER Genotype"))
 })
})
  


