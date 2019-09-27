library(shiny)
library(ggplot2)
#setwd("/Users/liz6/Box Sync/WCRC/Lee-Oesterreich-lab-personnel/Current/Vaciry/Vaciry/Bioinformatical analysis/Shiny/GEMM_Tumor/")
GEMM_CPM1=read.csv("GEMM_CPM1.csv",header=T, row.names=1)
GEMM_CPM3=read.csv("GEMM_CPM3.csv",header=T, row.names=1)
group<-read.csv("Group.csv",header=T,stringsAsFactors = F)
group2<-read.csv("Group2.csv",header=T,stringsAsFactors = F)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  output$Table<- renderTable(group)
  
  output$ER<- renderPlot({
    gene1<-input$gene
    GEMM_CPM1$ER<-group$ER
    df<-GEMM_CPM1[,c("ER",gene1)]
    colnames(df)<-c("ER","gene1")
    p<-wilcox.test(df[which(df$ER=="WT"),"gene1"], df[which(df$ER=="Mutant"),"gene1"], paired = F, alternative = "two.sided")$p.value
    p<-paste("Mann-Whitney U test=", 
             signif(p, digits=2))
    
    boxplot(gene1 ~ ER, data=df,outpch=NA,ylab="TMMlog2(CPM+1)",par(cex.lab=1.5),par(cex.axis=1.5), par(cex.main=2),xlab=p)
    stripchart(gene1 ~ ER, data = df, vertical = TRUE, method = "jitter",
               pch = 16, col = c("red","blue","black"),add = TRUE) 
    
    title(paste(gene1, "By Tumors"))
  })
  output$ER2<- renderPlot({
    gene1<-input$gene
    GEMM_CPM3$ER<-group2$ER_by_mouse
    df2<-GEMM_CPM3[,c("ER",gene1)]
    colnames(df2)<-c("ER","gene1")
    
    p2<-wilcox.test(df2[which(df2$ER=="WT"),"gene1"], df2[which(df2$ER=="Mutant"),"gene1"], paired = F, alternative = "two.sided")$p.value
    p2<-paste("Mann-Whitney U test=", 
              signif(p2, digits=2))
    
    boxplot(gene1 ~ ER, data=df2,outpch=NA,ylab="TMMlog2(CPM+1)", par(cex.lab=1.5),par(cex.axis=1.5),par(cex.main=2),xlab=p2)
    stripchart(gene1 ~ ER, data = df2, vertical = TRUE, method = "jitter",
               pch = 16, col = c("red","blue","black"),add = TRUE) 
    title(paste(gene1, "By Mouse"))
  })
})