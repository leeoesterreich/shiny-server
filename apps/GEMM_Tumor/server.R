library(shiny)
library(ggplot2)

GEMM_CPM<-read.csv("GEMM_CPM.csv")
group<-read.csv("Group.csv",header=T)
GEMM_CPM2<-read.csv("GEMM_CPM2.csv")
group2<-read.csv("Group2.csv",header=T)

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
   title(paste(gene1, "Expression~ER Genotype by Tumor"))
 })
 output$ER_by_mouse<- renderPlot({
   gene1<-input$gene
   rownames(GEMM_CPM2)=make.unique(as.character(GEMM_CPM2$Human))
   GEMM_CPM2=GEMM_CPM2[,2:15]
   GEMM_CPM2<-t(GEMM_CPM2)
   class(GEMM_CPM2)
   GEMM_CPM2=as.data.frame(GEMM_CPM2)
   GEMM_CPM2$ER_by_mouse<-group2$ER_by_mouse
   df2<-GEMM_CPM2[,c("ER_by_mouse",gene1)]
   colnames(df2)<-c("ER_by_mouse","gene1")
   
   p1<-wilcox.test(df2[which(df2$ER_by_mouse=="WT"),"gene1"], df2[which(df2$ER_by_mouse=="Mutant"),"gene1"], paired = F, alternative = "two.sided")$p.value
   p1<-paste("Mann-Whitney U p-value (WT vs Mut)=", signif(p2, digits=2))
   
   boxplot(gene1 ~ ER_by_mouse, data=df,outpch=NA,ylab="TMMLog2(CPM+1)", par(cex.lab=1.5),par(cex.axis=1.5),par(cex.main=2),xlab=p2)
   stripchart(gene1 ~ ER, data = df, vertical = TRUE, method = "jitter",
              pch = 20, col = c("red","black"),add = TRUE) 
   title(paste(gene1, "Expression~ER Genotype by Mouse"))
 })
})
  


