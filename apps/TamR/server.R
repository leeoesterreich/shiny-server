library(shiny)
library(ggplot2)
#setwd("/Users/vaciry/Box/WCRC/Lee-Oesterreich-lab-personnel/Current/Vaciry_Li/Vaciry/Bioinformatical analysis/Shiny/TamR/")
DF1<-read.csv("DF1.csv",sep=",",header=T)
DF2<-read.csv("DF2.csv",sep=",",header=T)
DF8<-read.csv("DF8.csv",sep=",",header=T)
DF4<-read.csv("DF4.csv",sep=",",header=T)
DF5<-read.csv("DF5.csv",sep=",",header=T)
DF6<-read.csv("DF6.csv",sep=",",header=T)
DF7<-read.csv("DF7.csv",sep=",",header=T)
DF11<-read.csv("DF11.csv",sep=",",header=T)
DF12<-read.csv("DF12.csv",sep=",",header=T)
Key<-read.csv("Key.csv",header=T)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  output$Table<- renderTable(Key)
  
  output$DF1 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), DF1[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("Resistance", "meanTPM", "sdTPM")
    rownames(df) = c("Parental","TamR")
    df = as.data.frame(df)
    
    df["Parental", "Resistance"] = "Parental"
    df["Parental", "meanTPM"] = mean(as.vector(as.matrix(DF1[row, 2:5])))
    df["Parental", "sdTPM"] = sd(as.vector(as.matrix(DF1[row, 2:5])))
    
    df["TamR", "Resistance"] = "TamR"
    df["TamR", "meanTPM"] = mean(as.vector(as.matrix(DF1[row, 6:9])))
    df["TamR", "sdTPM"] = sd(as.vector(as.matrix(DF1[row, 6:9])))
    
    df$Resistance = factor(df$Resistance, levels = c("Parental", "TamR"))
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = Resistance, y = meanTPM,
                               fill = factor(Resistance)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Model", y = "log2 (CPM+1)") +
      ggtitle(paste("MCF7_TamR_1")) +
      scale_fill_discrete(name = "Model") + 
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"), plot.title = element_text(size = 14),legend.position="none")
      
  })
  
  output$DF2 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), DF2[,"gene_id"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("Resistance", "meanTPM", "sdTPM")
    rownames(df) = c("Parental","TamR")
    df = as.data.frame(df)
    
    df["Parental", "Resistance"] = "Parental"
    df["Parental", "meanTPM"] = mean(as.vector(as.matrix(DF2[row, 2:4])))
    df["Parental", "sdTPM"] = sd(as.vector(as.matrix(DF2[row, 2:4])))
    
    df["TamR", "Resistance"] = "TamR"
    df["TamR", "meanTPM"] = mean(as.vector(as.matrix(DF2[row, 5:7])))
    df["TamR", "sdTPM"] = sd(as.vector(as.matrix(DF2[row, 5:7])))
    
    df$Resistance = factor(df$Resistance, levels = c("Parental", "TamR"))
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = Resistance, y = meanTPM,
                               fill = factor(Resistance)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Model", y = "log2 (CPM+1)") +
      ggtitle(paste("MCF7_TamR_2")) +
      scale_fill_discrete(name = "Model") + 
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"), plot.title = element_text(size = 14),legend.position="none")
    
  })
  
  output$DF4 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), DF4[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("Resistance", "meanTPM", "sdTPM")
    rownames(df) = c("Parental","TamR")
    df = as.data.frame(df)
    
    df["Parental", "Resistance"] = "Parental"
    df["Parental", "meanTPM"] = mean(as.vector(as.matrix(DF4[row, 2:4])))
    df["Parental", "sdTPM"] = sd(as.vector(as.matrix(DF4[row, 2:4])))
    
    df["TamR", "Resistance"] = "TamR"
    df["TamR", "meanTPM"] = mean(as.vector(as.matrix(DF4[row, 5:7])))
    df["TamR", "sdTPM"] = sd(as.vector(as.matrix(DF4[row, 5:7])))
    
    df$Resistance = factor(df$Resistance, levels = c("Parental", "TamR"))
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = Resistance, y = meanTPM,
                               fill = factor(Resistance)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Model", y = "log2 Intensity") +
      ggtitle(paste("MCF7_TamR_3")) +
      scale_fill_discrete(name = "Model") + 
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"), plot.title = element_text(size = 14),legend.position="none")
    
  })
  
  output$DF5 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), DF5[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("Resistance", "meanTPM", "sdTPM")
    rownames(df) = c("Parental","TamR")
    df = as.data.frame(df)
    
    df["Parental", "Resistance"] = "Parental"
    df["Parental", "meanTPM"] = mean(as.vector(as.matrix(DF5[row, 2:3])))
    df["Parental", "sdTPM"] = sd(as.vector(as.matrix(DF5[row, 2:3])))
    
    df["TamR", "Resistance"] = "TamR"
    df["TamR", "meanTPM"] = mean(as.vector(as.matrix(DF5[row, 6:7])))
    df["TamR", "sdTPM"] = sd(as.vector(as.matrix(DF5[row, 6:7])))
    
    df$Resistance = factor(df$Resistance, levels = c("Parental", "TamR"))
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = Resistance, y = meanTPM,
                               fill = factor(Resistance)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Model", y = "log2 (CPM+1)") +
      ggtitle(paste("MCF7_TamR_4")) +
      scale_fill_discrete(name = "Model") + 
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"), plot.title = element_text(size = 14),legend.position="none")
    
  })
  
  output$DF6 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), DF6[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("Resistance", "meanTPM", "sdTPM")
    rownames(df) = c("Parental","TamR")
    df = as.data.frame(df)
    
    df["Parental", "Resistance"] = "Parental"
    df["Parental", "meanTPM"] = mean(as.vector(as.matrix(DF6[row, 2:4])))
    df["Parental", "sdTPM"] = sd(as.vector(as.matrix(DF6[row, 2:4])))
    
    df["TamR", "Resistance"] = "TamR"
    df["TamR", "meanTPM"] = mean(as.vector(as.matrix(DF6[row, 5:7])))
    df["TamR", "sdTPM"] = sd(as.vector(as.matrix(DF6[row, 5:7])))
    
    df$Resistance = factor(df$Resistance, levels = c("Parental", "TamR"))
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = Resistance, y = meanTPM,
                               fill = factor(Resistance)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Model", y = "log2 (TPM+1)") +
      ggtitle(paste("MCF7_TamR_5")) +
      scale_fill_discrete(name = "Model") + 
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"), plot.title = element_text(size = 14),legend.position="none")
    
  })
  
  output$DF7_1 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), DF7[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("Resistance", "meanTPM", "sdTPM")
    rownames(df) = c("Parental","TamR")
    df = as.data.frame(df)
    
    df["Parental", "Resistance"] = "Parental"
    df["Parental", "meanTPM"] = mean(as.vector(as.matrix(DF7[row, 2])))
    df["Parental", "sdTPM"] = sd(as.vector(as.matrix(DF7[row, 2])))
    
    df["TamR", "Resistance"] = "TamR"
    df["TamR", "meanTPM"] = mean(as.vector(as.matrix(DF7[row, 3])))
    df["TamR", "sdTPM"] = sd(as.vector(as.matrix(DF7[row, 3])))
    
    df$Resistance = factor(df$Resistance, levels = c("Parental", "TamR"))
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = Resistance, y = meanTPM,
                               fill = factor(Resistance)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Model", y = "log2 (CPM+1)") +
      ggtitle(paste("MCF7_TamR_6")) +
      scale_fill_discrete(name = "Model") + 
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"), plot.title = element_text(size = 14),legend.position="none")
    
  })
 
  output$DF7_2 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), DF7[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("Resistance", "meanTPM", "sdTPM")
    rownames(df) = c("Parental","TamR")
    df = as.data.frame(df)
    
    df["Parental", "Resistance"] = "Parental"
    df["Parental", "meanTPM"] = mean(as.vector(as.matrix(DF7[row, 4])))
    df["Parental", "sdTPM"] = sd(as.vector(as.matrix(DF7[row, 4])))
    
    df["TamR", "Resistance"] = "TamR"
    df["TamR", "meanTPM"] = mean(as.vector(as.matrix(DF7[row, 5:6])))
    df["TamR", "sdTPM"] = sd(as.vector(as.matrix(DF7[row, 5:6])))
    
    df$Resistance = factor(df$Resistance, levels = c("Parental", "TamR"))
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = Resistance, y = meanTPM,
                               fill = factor(Resistance)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Model", y = "log2 (CPM+1)") +
      ggtitle(paste("T47D_TamR")) +
      scale_fill_discrete(name = "Model") + 
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"), plot.title = element_text(size = 14),legend.position="none")
    
  })
  output$DF7_3 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), DF7[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("Resistance", "meanTPM", "sdTPM")
    rownames(df) = c("Parental","TamR")
    df = as.data.frame(df)
    
    df["Parental", "Resistance"] = "Parental"
    df["Parental", "meanTPM"] = mean(as.vector(as.matrix(DF7[row, 7])))
    df["Parental", "sdTPM"] = sd(as.vector(as.matrix(DF7[row, 7])))
    
    df["TamR", "Resistance"] = "TamR"
    df["TamR", "meanTPM"] = mean(as.vector(as.matrix(DF7[row, 8:9])))
    df["TamR", "sdTPM"] = sd(as.vector(as.matrix(DF7[row, 8:9])))
    
    df$Resistance = factor(df$Resistance, levels = c("Parental", "TamR"))
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = Resistance, y = meanTPM,
                               fill = factor(Resistance)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Model", y = "log2 (CPM+1)") +
      ggtitle(paste("ZR75-1_TamR")) +
      scale_fill_discrete(name = "Model") + 
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"), plot.title = element_text(size = 14),legend.position="none")
    
  })
  output$DF7_4 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), DF7[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("Resistance", "meanTPM", "sdTPM")
    rownames(df) = c("Parental","TamR")
    df = as.data.frame(df)
    
    df["Parental", "Resistance"] = "Parental"
    df["Parental", "meanTPM"] = mean(as.vector(as.matrix(DF7[row, 10])))
    df["Parental", "sdTPM"] = sd(as.vector(as.matrix(DF7[row, 10])))
    
    df["TamR", "Resistance"] = "TamR"
    df["TamR", "meanTPM"] = mean(as.vector(as.matrix(DF7[row, 11:12])))
    df["TamR", "sdTPM"] = sd(as.vector(as.matrix(DF7[row, 10:12])))
    
    df$Resistance = factor(df$Resistance, levels = c("Parental", "TamR"))
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = Resistance, y = meanTPM,
                               fill = factor(Resistance)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Model", y = "log2 (CPM+1)") +
      ggtitle(paste("BT474_TamR")) +
      scale_fill_discrete(name = "Model") + 
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"), plot.title = element_text(size = 14),legend.position="none")
    
  })
  output$DF8 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), DF8[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("Resistance", "meanTPM", "sdTPM")
    rownames(df) = c("Parental","TamR")
    df = as.data.frame(df)
    
    df["Parental", "Resistance"] = "Parental"
    df["Parental", "meanTPM"] = mean(as.vector(as.matrix(DF8[row, 2:4])))
    df["Parental", "sdTPM"] = sd(as.vector(as.matrix(DF8[row, 2:4])))
    
    df["TamR", "Resistance"] = "TamR"
    df["TamR", "meanTPM"] = mean(as.vector(as.matrix(DF8[row, 5:7])))
    df["TamR", "sdTPM"] = sd(as.vector(as.matrix(DF8[row, 5:7])))
    
    df$Resistance = factor(df$Resistance, levels = c("Parental", "TamR"))
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = Resistance, y = meanTPM,
                               fill = factor(Resistance)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Model", y = "log2 Intensity") +
      ggtitle(paste("SUM44_TamR")) +
      scale_fill_discrete(name = "Model") + 
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"), plot.title = element_text(size = 14),legend.position="none")
    
  })
  output$DF11 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), DF11[,"gene_id"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("Resistance", "meanTPM", "sdTPM")
    rownames(df) = c("Parental","FasR")
    df = as.data.frame(df)
    
    df["Parental", "Resistance"] = "Parental"
    df["Parental", "meanTPM"] = mean(as.vector(as.matrix(DF11[row, 2:4])))
    df["Parental", "sdTPM"] = sd(as.vector(as.matrix(DF11[row, 2:4])))
    
    df["FasR", "Resistance"] = "FasR"
    df["FasR", "meanTPM"] = mean(as.vector(as.matrix(DF11[row, 5:7])))
    df["FasR", "sdTPM"] = sd(as.vector(as.matrix(DF11[row, 5:7])))
    
    df$Resistance = factor(df$Resistance, levels = c("Parental", "FasR"))
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = Resistance, y = meanTPM,
                               fill = factor(Resistance)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Model", y = "log2 (CPM+1)") +
      ggtitle(paste("MCF7_FasR_1")) +
      scale_fill_discrete(name = "Model") + 
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"), plot.title = element_text(size = 14),legend.position="none")
    
  })
  output$DF12 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), DF12[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("Resistance", "meanTPM", "sdTPM")
    rownames(df) = c("Parental","FasR")
    df = as.data.frame(df)
    
    df["Parental", "Resistance"] = "Parental"
    df["Parental", "meanTPM"] = mean(as.vector(as.matrix(DF12[row, 2:4])))
    df["Parental", "sdTPM"] = sd(as.vector(as.matrix(DF12[row, 2:4])))
    
    df["FasR", "Resistance"] = "FasR"
    df["FasR", "meanTPM"] = mean(as.vector(as.matrix(DF12[row, 5:7])))
    df["FasR", "sdTPM"] = sd(as.vector(as.matrix(DF12[row, 5:7])))
    
    df$Resistance = factor(df$Resistance, levels = c("Parental", "FasR"))
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = Resistance, y = meanTPM,
                               fill = factor(Resistance)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Model", y = "log2 (TPM+1)") +
      ggtitle(paste("MCF7_FasR_2")) +
      scale_fill_discrete(name = "Model") + 
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14,face="bold"), plot.title = element_text(size = 14),legend.position="none")
    
  })
})
