library(shiny)
library(ggplot2)

MCF7_1<-read.csv("SO_MCF7.csv",sep=",",header=T)
MCF7_2<-read.csv("JG_MCF7.csv",sep=",",header=T)
MCF7_3<-read.csv("SA_MCF7.csv",sep=",",header=T)
MCF7_4<-read.csv("MB_MCF7.csv",sep=",",header=T)
T47D_1<-read.csv("SO_T47D.csv",sep=",",header=T)
T47D_2<-read.csv("Creighton.csv",sep=",",header=T)
T47D_3<-read.csv("JG_T47D.csv",sep=",",header=T)
T47D_4<-read.csv("Shapiro_T47D.csv",sep=",",header=T)
ZR<-read.csv("ZR75_1_Log2INT.csv",sep=",",header=T)
Sikora<-read.csv("Sikora.csv",sep=",",header=T)
Key<-read.csv("Key.csv",header=T)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  output$Table<- renderTable(Key)
  
  output$MCF7_1 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), MCF7_1[,"Hugo.name"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("E2_Veh", "meanTPM", "sdTPM")
    rownames(df) = c("Veh","E2")
    df = as.data.frame(df)
    
    df["Veh", "E2_Veh"] = "Veh"
    df["Veh", "meanTPM"] = mean(as.vector(as.matrix(MCF7_1[row, 2:5])))
    df["Veh", "sdTPM"] = sd(as.vector(as.matrix(MCF7_1[row, 2:5])))
    
    df["E2", "E2_Veh"] = "E2"
    df["E2", "meanTPM"] = mean(as.vector(as.matrix(MCF7_1[row, 6:9])))
    df["E2", "sdTPM"] = sd(as.vector(as.matrix(MCF7_1[row, 6:9])))
    
    
    df$E2_Veh = factor(df$E2_Veh, levels = c("Veh", "E2"))
    
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = E2_Veh, y = meanTPM,
                               fill = factor(E2_Veh)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Treatment", y = "log2 (TPM+1)") +
      ggtitle(paste("MCF7_1", input$gene)) +
      scale_fill_discrete(name = "Treatments") + 
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"), plot.title = element_text(size = rel(2)),legend.position="none")
      
  })
  
  output$MCF7_2 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), MCF7_2[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("E2_Veh", "meanTPM", "sdTPM")
    rownames(df) = c("Veh","E2")
    df = as.data.frame(df)
    
    df["Veh", "E2_Veh"] = "Veh"
    df["Veh", "meanTPM"] = mean(as.vector(as.matrix(MCF7_2[row, 2:3])))
    df["Veh", "sdTPM"] = sd(as.vector(as.matrix(MCF7_2[row, 2:3])))
    
    df["E2", "E2_Veh"] = "E2"
    df["E2", "meanTPM"] = mean(as.vector(as.matrix(MCF7_2[row, 4:5])))
    df["E2", "sdTPM"] = sd(as.vector(as.matrix(MCF7_2[row, 4:5])))
    
    
    df$E2_Veh = factor(df$E2_Veh, levels = c("Veh", "E2"))
    
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = E2_Veh, y = meanTPM,
                               fill = factor(E2_Veh)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Treatment", y = "log2 (CPM+1)") +
      ggtitle(paste("MCF7_2", input$gene)) +
      scale_fill_discrete(name = "Treatments") + 
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"), plot.title = element_text(size = rel(2)),legend.position="none")
    
  })
  
  output$MCF7_3 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), MCF7_3[,"hugo"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("E2_Veh", "meanTPM", "sdTPM")
    rownames(df) = c("Veh","E2")
    df = as.data.frame(df)
    
    df["Veh", "E2_Veh"] = "Veh"
    df["Veh", "meanTPM"] = mean(as.vector(as.matrix(MCF7_3[row, 2:4])))
    df["Veh", "sdTPM"] = sd(as.vector(as.matrix(MCF7_3[row, 2:4])))
    
    df["E2", "E2_Veh"] = "E2"
    df["E2", "meanTPM"] = mean(as.vector(as.matrix(MCF7_3[row, 5:7])))
    df["E2", "sdTPM"] = sd(as.vector(as.matrix(MCF7_3[row, 5:7])))
    
    
    df$E2_Veh = factor(df$E2_Veh, levels = c("Veh", "E2"))
    
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = E2_Veh, y = meanTPM,
                               fill = factor(E2_Veh)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Treatment", y = "log2 (TPM+1)") +
      ggtitle(paste("MCF7_3", input$gene)) +
      scale_fill_discrete(name = "Treatments") + 
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"), plot.title = element_text(size = rel(2)),legend.position="none")
    
  })
  output$MCF7_4 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), MCF7_4[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("E2_Veh", "meanTPM", "sdTPM")
    rownames(df) = c("Veh","E2")
    df = as.data.frame(df)
    
    df["Veh", "E2_Veh"] = "Veh"
    df["Veh", "meanTPM"] = mean(as.vector(as.matrix(MCF7_4[row, 2:4])))
    df["Veh", "sdTPM"] = sd(as.vector(as.matrix(MCF7_4[row, 2:4])))
    
    df["E2", "E2_Veh"] = "E2"
    df["E2", "meanTPM"] = mean(as.vector(as.matrix(MCF7_4[row, 5:7])))
    df["E2", "sdTPM"] = sd(as.vector(as.matrix(MCF7_4[row, 5:7])))
    
    
    df$E2_Veh = factor(df$E2_Veh, levels = c("Veh", "E2"))
    
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = E2_Veh, y = meanTPM,
                               fill = factor(E2_Veh)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Treatment", y = "log2 (CPM+1)") +
      ggtitle(paste("MCF7_4", input$gene)) +
      scale_fill_discrete(name = "Treatments") + 
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"), plot.title = element_text(size = rel(2)),legend.position="none")
    
  })
  
  output$T47D_1 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), T47D_1[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("E2_Veh", "meanTPM", "sdTPM")
    rownames(df) = c("Veh","E2")
    df = as.data.frame(df)
    
    df["Veh", "E2_Veh"] = "Veh"
    df["Veh", "meanTPM"] = mean(as.vector(as.matrix(T47D_1[row, 2:5])))
    df["Veh", "sdTPM"] = sd(as.vector(as.matrix(T47D_1[row, 2:5])))
    
    df["E2", "E2_Veh"] = "E2"
    df["E2", "meanTPM"] = mean(as.vector(as.matrix(T47D_1[row, 6:9])))
    df["E2", "sdTPM"] = sd(as.vector(as.matrix(T47D_1[row, 6:9])))
    
    
    df$E2_Veh = factor(df$E2_Veh, levels = c("Veh", "E2"))
    
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = E2_Veh, y = meanTPM,
                               fill = factor(E2_Veh)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Treatment", y = "Log2 (CPM+1)") +
      ggtitle(paste("T47D_1", input$gene)) +
      scale_fill_discrete(name = "Treatments") + 
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"), plot.title = element_text(size = rel(2)),legend.position="none")
    
  })
  output$T47D_2 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), T47D_2[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("E2_Veh", "meanTPM", "sdTPM")
    rownames(df) = c("Veh","E2")
    df = as.data.frame(df)
    
    df["Veh", "E2_Veh"] = "Veh"
    df["Veh", "meanTPM"] = mean(as.vector(as.matrix(T47D_2[row, 2:3])))
    df["Veh", "sdTPM"] = sd(as.vector(as.matrix(T47D_2[row, 2:3])))
    
    df["E2", "E2_Veh"] = "E2"
    df["E2", "meanTPM"] = mean(as.vector(as.matrix(T47D_2[row, 4:5])))
    df["E2", "sdTPM"] = sd(as.vector(as.matrix(T47D_2[row, 4:5])))
    
    
    df$E2_Veh = factor(df$E2_Veh, levels = c("Veh", "E2"))
    
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = E2_Veh, y = meanTPM,
                               fill = factor(E2_Veh)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Treatment", y = "Log2 Norm Intensity") +
      ggtitle(paste("T47D_2", input$gene)) +
      scale_fill_discrete(name = "Treatments") + 
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"), plot.title = element_text(size = rel(2)),legend.position="none")
    
  })
 
  output$T47D_3 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), T47D_3[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("E2_Veh", "meanTPM", "sdTPM")
    rownames(df) = c("Veh","E2")
    df = as.data.frame(df)
    
    df["Veh", "E2_Veh"] = "Veh"
    df["Veh", "meanTPM"] = mean(as.vector(as.matrix(T47D_3[row, 2:3])))
    df["Veh", "sdTPM"] = sd(as.vector(as.matrix(T47D_3[row, 2:3])))
    
    df["E2", "E2_Veh"] = "E2"
    df["E2", "meanTPM"] = mean(as.vector(as.matrix(T47D_3[row, 4:5])))
    df["E2", "sdTPM"] = sd(as.vector(as.matrix(T47D_3[row, 4:5])))
    
    
    df$E2_Veh = factor(df$E2_Veh, levels = c("Veh", "E2"))
    
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = E2_Veh, y = meanTPM,
                               fill = factor(E2_Veh)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Treatment", y = "Log2 (CPM+1)") +
      ggtitle(paste("T47D_3", input$gene)) +
      scale_fill_discrete(name = "Treatments") + 
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"), plot.title = element_text(size = rel(2)),legend.position="none")
    
  })
  output$T47D_4 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), T47D_4[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("E2_Veh", "meanTPM", "sdTPM")
    rownames(df) = c("Veh","E2")
    df = as.data.frame(df)
    
    df["Veh", "E2_Veh"] = "Veh"
    df["Veh", "meanTPM"] = mean(as.vector(as.matrix(T47D_4[row, 2:4])))
    df["Veh", "sdTPM"] = sd(as.vector(as.matrix(T47D_4[row, 2:4])))
    
    df["E2", "E2_Veh"] = "E2"
    df["E2", "meanTPM"] = mean(as.vector(as.matrix(T47D_4[row, 5:7])))
    df["E2", "sdTPM"] = sd(as.vector(as.matrix(T47D_4[row, 5:7])))
    
    
    df$E2_Veh = factor(df$E2_Veh, levels = c("Veh", "E2"))
    
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = E2_Veh, y = meanTPM,
                               fill = factor(E2_Veh)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Treatment", y = "Log2 (CPM+1)") +
      ggtitle(paste("T47D_4", input$gene)) +
      scale_fill_discrete(name = "Treatments") + 
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"), plot.title = element_text(size = rel(2)),legend.position="none")
    
  })
  output$BT474 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), T47D_2[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("E2_Veh", "meanTPM", "sdTPM")
    rownames(df) = c("Veh","E2")
    df = as.data.frame(df)
    
    df["Veh", "E2_Veh"] = "Veh"
    df["Veh", "meanTPM"] = mean(as.vector(as.matrix(T47D_2[row, 6:7])))
    df["Veh", "sdTPM"] = sd(as.vector(as.matrix(T47D_2[row, 6:7])))
    
    df["E2", "E2_Veh"] = "E2"
    df["E2", "meanTPM"] = mean(as.vector(as.matrix(T47D_2[row, 8:9])))
    df["E2", "sdTPM"] = sd(as.vector(as.matrix(T47D_2[row, 8:9])))
    
    
    df$E2_Veh = factor(df$E2_Veh, levels = c("Veh", "E2"))
    
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = E2_Veh, y = meanTPM,
                               fill = factor(E2_Veh)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Treatment", y = "Log2 Norm Intensity") +
      ggtitle(paste("BT474", input$gene)) +
      scale_fill_discrete(name = "Treatments") + 
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"), plot.title = element_text(size = rel(2)),legend.position="none")
    
  })
  output$ZR <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), ZR[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("E2_Veh", "meanTPM", "sdTPM")
    rownames(df) = c("Veh","E2")
    df = as.data.frame(df)
    
    df["Veh", "E2_Veh"] = "Veh"
    df["Veh", "meanTPM"] = mean(as.vector(as.matrix(ZR[row, 2:5])))
    df["Veh", "sdTPM"] = sd(as.vector(as.matrix(ZR[row, 2:5])))
    
    df["E2", "E2_Veh"] = "E2"
    df["E2", "meanTPM"] = mean(as.vector(as.matrix(ZR[row, 6:9])))
    df["E2", "sdTPM"] = sd(as.vector(as.matrix(ZR[row, 6:9])))
    
    
    df$E2_Veh = factor(df$E2_Veh, levels = c("Veh", "E2"))
    
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = E2_Veh, y = meanTPM,
                               fill = factor(E2_Veh)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Treatment", y = "Log2 Norm Intensity") +
      ggtitle(paste("ZR75-1", input$gene)) +
      scale_fill_discrete(name = "Treatments") + 
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"), plot.title = element_text(size = rel(2)),legend.position="none")
    
  })
  output$MM134 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), Sikora[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("E2_Veh", "meanTPM", "sdTPM")
    rownames(df) = c("Veh","E2")
    df = as.data.frame(df)
    
    df["Veh", "E2_Veh"] = "Veh"
    df["Veh", "meanTPM"] = mean(as.vector(as.matrix(Sikora[row, 2:5])))
    df["Veh", "sdTPM"] = sd(as.vector(as.matrix(Sikora[row, 2:5])))
    
    df["E2", "E2_Veh"] = "E2"
    df["E2", "meanTPM"] = mean(as.vector(as.matrix(Sikora[row, 6:9])))
    df["E2", "sdTPM"] = sd(as.vector(as.matrix(Sikora[row, 6:9])))
    
    
    df$E2_Veh = factor(df$E2_Veh, levels = c("Veh", "E2"))
    
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = E2_Veh, y = meanTPM,
                               fill = factor(E2_Veh)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Treatment", y = "Log2 Norm Intensity") +
      ggtitle(paste("MM134", input$gene)) +
      scale_fill_discrete(name = "Treatments") + 
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"), plot.title = element_text(size = rel(2)),legend.position="none")
    
  })
  output$SUM44 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), Sikora[,"X"])
    df = matrix(NA, nrow = 2, ncol =  3)
    colnames(df) = c("E2_Veh", "meanTPM", "sdTPM")
    rownames(df) = c("Veh","E2")
    df = as.data.frame(df)
    
    df["Veh", "E2_Veh"] = "Veh"
    df["Veh", "meanTPM"] = mean(as.vector(as.matrix(Sikora[row, 10:13])))
    df["Veh", "sdTPM"] = sd(as.vector(as.matrix(Sikora[row, 10:13])))
    
    df["E2", "E2_Veh"] = "E2"
    df["E2", "meanTPM"] = mean(as.vector(as.matrix(Sikora[row, 14:17])))
    df["E2", "sdTPM"] = sd(as.vector(as.matrix(Sikora[row, 14:17])))
    
    
    df$E2_Veh = factor(df$E2_Veh, levels = c("Veh", "E2"))
    
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = E2_Veh, y = meanTPM,
                               fill = factor(E2_Veh)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9),fill=c("grey","red")) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "Treatment", y = "Log2 Norm Intensity") +
      ggtitle(paste("SUM44", input$gene)) +
      scale_fill_discrete(name = "Treatments") + 
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"), plot.title = element_text(size = rel(2)),legend.position="none")
    
  })
})
