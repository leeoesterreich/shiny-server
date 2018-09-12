library(shiny)
library(ggplot2)


#load("esr1_salmonlog2TPM_HGNC.Rda")
#esr1_salmonlog2TPM_HGNC = temp

load("esr1_salmonlog2TPM_HGNC.Rda")
esr1_salmonlog2TPM_HGNC =tmp

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  output$plotTPM_MCF7 <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), esr1_salmonlog2TPM_HGNC[,"external_gene_name"])
    df = matrix(NA, nrow = 6, ncol =  4)
    colnames(df) = c("WT_MUT", "E2_Veh", "meanTPM", "sdTPM")
    rownames(df) = c("WT_Veh","WT_E2","Y537S_Veh","Y537S_E2","D538G_Veh", "D538G_E2")
    df = as.data.frame(df)
    
    df["WT_Veh", "WT_MUT"] = "WT"
    df["WT_Veh", "E2_Veh"] = "Veh"
    df["WT_Veh", "meanTPM"] = mean(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 2:5])))
    df["WT_Veh", "sdTPM"] = sd(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 2:5])))
    
    df["WT_E2", "WT_MUT"] = "WT"
    df["WT_E2", "E2_Veh"] = "E2"
    df["WT_E2", "meanTPM"] = mean(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 14:17])))
    df["WT_E2", "sdTPM"] = sd(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 14:17])))
    
    df["Y537S_Veh", "WT_MUT"] = "Y537S"
    df["Y537S_Veh", "E2_Veh"] = "Veh"
    df["Y537S_Veh", "meanTPM"] = mean(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 6:9])))
    df["Y537S_Veh", "sdTPM"] = sd(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 6:9])))
    
    df["Y537S_E2", "WT_MUT"] = "Y537S"
    df["Y537S_E2", "E2_Veh"] = "E2"
    df["Y537S_E2", "meanTPM"] = mean(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 18:21])))
    df["Y537S_E2", "sdTPM"] = sd(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 18:21])))
    
    df["D538G_Veh", "WT_MUT"] = "D538G"
    df["D538G_Veh", "E2_Veh"] = "Veh"
    df["D538G_Veh", "meanTPM"] = mean(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 10:13])))
    df["D538G_Veh", "sdTPM"] = sd(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 10:13])))
    
    df["D538G_E2", "WT_MUT"] = "D538G"
    df["D538G_E2", "E2_Veh"] = "E2"
    df["D538G_E2", "meanTPM"] = mean(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 22:25])))
    df["D538G_E2", "sdTPM"] = sd(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 22:25])))
    
    df$WT_MUT = factor(df$WT_MUT, levels = c("WT", "Y537S", "D538G"))
    df$E2_Veh = factor(df$E2_Veh, levels = c("Veh", "E2"))
    
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = WT_MUT, y = meanTPM,
                                   fill = factor(E2_Veh)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9)) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "ER Genotype", y = "log2 TPM") +
      ggtitle(paste("MCF7 log2 expression of", input$gene)) +
      scale_fill_discrete(name = "E2 status") + theme(axis.text=element_text(size=18),
           axis.title=element_text(size=20,face="bold"), plot.title = element_text(size = rel(2)), legend.text = element_text(size = 18), legend.title = element_text(size=20))
      #scale_fill_discrete(name = "E2 status") + theme(axis.text=element_text(size=9),
      #                                                axis.title=element_text(size=10,face="bold"), plot.title = element_text(size = 10), legend.text = element_text(size = 9), legend.title = element_text(size=10))
      })
  
  output$plotTPM_T47D <- renderPlot({
    row = grep(paste("^", input$gene, "$", sep=""), esr1_salmonlog2TPM_HGNC[,"external_gene_name"])
    df = matrix(NA, nrow = 6, ncol =  4)
    colnames(df) = c("WT_MUT", "E2_Veh", "meanTPM", "sdTPM")
    rownames(df) = c("WT_Veh","WT_E2","Y537S_Veh","Y537S_E2","D538G_Veh", "D538G_E2")
    df = as.data.frame(df)
    
    df["WT_Veh", "WT_MUT"] = "WT"
    df["WT_Veh", "E2_Veh"] = "Veh"
    df["WT_Veh", "meanTPM"] = mean(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 26:29])))
    df["WT_Veh", "sdTPM"] = sd(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 26:29])))
    
    df["WT_E2", "WT_MUT"] = "WT"
    df["WT_E2", "E2_Veh"] = "E2"
    df["WT_E2", "meanTPM"] = mean(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 38:41])))
    df["WT_E2", "sdTPM"] = sd(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 38:41])))
    
    df["Y537S_Veh", "WT_MUT"] = "Y537S"
    df["Y537S_Veh", "E2_Veh"] = "Veh"
    df["Y537S_Veh", "meanTPM"] = mean(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 30:33])))
    df["Y537S_Veh", "sdTPM"] = sd(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 30:33])))
    
    df["Y537S_E2", "WT_MUT"] = "Y537S"
    df["Y537S_E2", "E2_Veh"] = "E2"
    df["Y537S_E2", "meanTPM"] = mean(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 42:45])))
    df["Y537S_E2", "sdTPM"] = sd(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 42:45])))
    
    df["D538G_Veh", "WT_MUT"] = "D538G"
    df["D538G_Veh", "E2_Veh"] = "Veh"
    df["D538G_Veh", "meanTPM"] = mean(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 34:37])))
    df["D538G_Veh", "sdTPM"] = sd(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 34:37])))
    
    df["D538G_E2", "WT_MUT"] = "D538G"
    df["D538G_E2", "E2_Veh"] = "E2"
    df["D538G_E2", "meanTPM"] = mean(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 46:49])))
    df["D538G_E2", "sdTPM"] = sd(as.vector(as.matrix(esr1_salmonlog2TPM_HGNC[row, 46:49])))
    
    df$WT_MUT = factor(df$WT_MUT, levels = c("WT", "Y537S", "D538G"))
    df$E2_Veh = factor(df$E2_Veh, levels = c("Veh", "E2"))
    
    
    #dodge = spacing between bars 
    limits <- aes(ymax = df$meanTPM + df$sdTPM,
                  ymin = df$meanTPM - df$sdTPM)
    
    p <- ggplot(data = df, aes(x = WT_MUT, y = meanTPM,
                               fill = factor(E2_Veh)))
    
    p + geom_bar(stat = "identity",
                 position = position_dodge(0.9)) +
      geom_errorbar(limits, position = position_dodge(0.9),
                    width = 0.25) +
      labs(x = "ER Genotype", y = "log2 TPM") +
      ggtitle(paste("T47D log2 expression of", input$gene)) +
      scale_fill_discrete(name = "E2 status") + theme(axis.text=element_text(size=18),
                                                      axis.title=element_text(size=20,face="bold"), plot.title = element_text(size = rel(2)), legend.text = element_text(size = 18), legend.title = element_text(size=20))
  })
  
 
})
