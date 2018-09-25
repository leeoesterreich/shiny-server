library(shiny)
library(ggplot2)
#load(PDX_counts)

HCI_CPM<-read.csv("PDX_CPM.csv")
group<-read.csv("Group.csv",stringsAsFactors = F)


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  output$Table<- renderTable(group)
   
  output$Types<- renderPlot({
   gene1<-input$gene
   rownames(HCI_CPM)=HCI_CPM$X
   HCI_CPM=HCI_CPM[,2:23]
   HCI_CPM1<-t(HCI_CPM)
   class(HCI_CPM1)
   HCI_CPM1=as.data.frame(HCI_CPM1)
   HCI_CPM1$Type<-group$Type
   df<-HCI_CPM1[,c("Type",gene1)]
   colnames(df)<-c("Type","gene1")
   
   p<-wilcox.test(df[which(df$Type=="ER+"),"gene1"], df[which(df$Type=="TNBC"),"gene1"], paired = F, alternative = "two.sided")$p.value
   p<-paste("Mann-Whitney U p-value (TNBC vs ER+)=", signif(p, digits=2))
  
   boxplot(gene1 ~ Type, data=df,outpch=NA,ylab="log2CPM",par(cex.lab=1.5),par(cex.axis=1.5), par(cex.main=2),xlab=p)
   stripchart(gene1 ~ Type, data = df, vertical = TRUE, method = "jitter",
              pch = 20, col = c("red","blue","black"),add = TRUE) 

   title(paste(gene1, "Expression_BrCA subtypes"))
  })
 output$ER<- renderPlot({
   gene1<-input$gene
   rownames(HCI_CPM)=HCI_CPM$X
   HCI_CPM=HCI_CPM[,2:23]
   HCI_CPM1<-t(HCI_CPM)
   class(HCI_CPM1)
   HCI_CPM1=as.data.frame(HCI_CPM1)
   HCI_CPM1$ER<-group$ER
   df<-HCI_CPM1[,c("ER",gene1)]
   colnames(df)<-c("ER","gene1")
   
   p2<-wilcox.test(df[which(df$ER=="WT ER"),"gene1"], df[which(df$ER=="Mut ER"),"gene1"], paired = F, alternative = "two.sided")$p.value
   p2<-paste("Mann-Whitney U p-value (WT vs Mut)=", signif(p2, digits=2))
   
   boxplot(gene1 ~ ER, data=df,outpch=NA,ylab="log2CPM", par(cex.lab=1.5),par(cex.axis=1.5),par(cex.main=2),xlab=p2)
   stripchart(gene1 ~ ER, data = df, vertical = TRUE, method = "jitter",
              pch = 20, col = c("red","blue","black"),add = TRUE) 
   title(paste(gene1, "Expression_ER Status"))
 })
})
  


