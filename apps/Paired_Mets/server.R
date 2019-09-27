library(shiny)
library(ggplot2)
library(data.table)
library(reshape)
#setwd("C:/Users/basudana/Box Sync/R/RNAseq_Brian_Bone_Ovaries_Salmon0.8.2/CPM/")
#load("PanMets.log2.tmmNorm.cpm.salmon.0.8.2.Rda")

#keep samples + gene column only
#PanMets <- tmp[,c(2:136)] #60619   135
#setwd("C:/Users/basudana/Box Sync/R/RNAseq_Brian_Bone_Ovaries_Salmon0.8.2/shinyapp")
#save(PanMets, file = "ShinyApp.PanMets.log2.tmmNorm.cpm.salmon.0.8.2.Rda")

#setwd(setwd("~/Desktop/shinyapp2"))
load("ShinyApp.PanMets.log2.tmmNorm.cpm.salmon.0.8.2.Rda")

Key_file= read.csv("SampleAnnotation.csv")

shinyServer(function(input, output) {
  
  #ovaries only
  output$ov <-renderPlot({
    gene<-input$gene
    chosedGene <- subset(PanMets, (external_gene_name %in% gene))
    
    mdat = melt(chosedGene, id.vars=("external_gene_name"))
    
    mdat_wkey<- merge(mdat,Key_file,by.x = "variable",by.y = "ID")
    colnames(mdat_wkey)[3] <- "log2CPM"
    
    mdat_tumors <- mdat_wkey
    #PAIRS
    mdat_pairs <- mdat_tumors[!(is.na(mdat_tumors$pair) | mdat_tumors$pair==""),]
    mdat_pairs <- mdat_pairs[!(mdat_pairs$subtype=="unknown" |  mdat_pairs$pair==""),]
    mdat_pairs$group <- factor(mdat_pairs$group, c("Primary", "MET")) #change order so primary appear left on plot
    mdat_pairs$site <- factor(mdat_pairs$site, c("Brain","Bone", "Ovary","Breast","GI"))
    mdat_pairs$ER_status <- factor(mdat_pairs$ER_status, c("Pos","Neg","Unk"))
    ovaries <- mdat_pairs[mdat_pairs$code== "Ov",]
    
    test<- subset(ovaries, (external_gene_name %in% gene))
    test.primary <- subset(test, (group %in% "Primary"))
    test.met <- subset(test, (group %in% "MET"))
    w<-wilcox.test(test.primary$log2CPM, test.met$log2CPM,paired = TRUE, alternative = "two.sided")
    
    ggplot(subset(ovaries, (external_gene_name %in% gene)), aes(x=group, y=log2CPM, color= subtype,group=pair))  +
      geom_point(size=3) +geom_line(aes(linetype=ER_status),size=0.5) + theme(axis.text=element_text(size=15,face="bold"),
                                                                              axis.title=element_text(size=14,face="bold")) +
      ggtitle(paste(gene,"Expression in Paired Ovarian Metastases")) + 
      theme(plot.title = element_text(color="black", face="bold", size=14, hjust=.5,vjust=1)) +
      labs(x="",y="log2(TMM-normalized CPM +1)",caption = paste("Paired Wilcoxon.test p value=",format(w$p.value, scientific = TRUE))) +
      scale_color_manual(values= c("purple4", "red4","darkorange")) +
      scale_linetype_manual(values= c(1,3)) +
      theme(legend.title=element_text(size=14,face="bold"))+
      theme(legend.text=element_text(size=12,face="bold"))
  })
  
  output$ov.sb <-renderPlot({
    gene<-input$gene
    chosedGene <- subset(PanMets, (external_gene_name %in% gene))
    
    mdat = melt(chosedGene, id.vars=("external_gene_name"))
    
    mdat_wkey<- merge(mdat,Key_file,by.x = "variable",by.y = "ID")
    colnames(mdat_wkey)[3] <- "log2CPM"
    
    mdat_tumors <- mdat_wkey
    #PAIRS
    mdat_pairs <- mdat_tumors[!(is.na(mdat_tumors$pair) | mdat_tumors$pair==""),]
    mdat_pairs <- mdat_pairs[!(mdat_pairs$subtype=="unknown" |  mdat_pairs$pair==""),]
    mdat_pairs$group <- factor(mdat_pairs$group, c("Primary", "MET")) #change order so primary appear left on plot
    mdat_pairs$site <- factor(mdat_pairs$site, c("Brain","Bone", "Ovary","Breast","GI"))
    mdat_pairs$ER_status <- factor(mdat_pairs$ER_status, c("Pos","Neg","Unk"))
    ova.all <- subset(mdat_tumors, code %in% "Ov")
    ovaries.MET <- ova.all[ova.all$group== "MET",]
    ovaries.MET.gene <- subset(ovaries.MET, external_gene_name %in% gene)
    ovaries.MET.ILC <- ovaries.MET.gene[ovaries.MET.gene$subtype=="ILC",]
    ovaries.MET.IDC <- ovaries.MET.gene[ovaries.MET.gene$subtype=="IDC",]
    w2<-wilcox.test(ovaries.MET.IDC$log2CPM,ovaries.MET.ILC$log2CPM,paired = FALSE, alternative = "two.sided")
    
    boxplot(ovaries.MET.IDC$log2CPM,ovaries.MET.ILC$log2CPM,main= paste(gene,"expression in Ovarian Mets(IDC vs ILC)"),cex.main=1.5,outline = TRUE,pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5),
            log = "",
            xlab=paste("Non-Paired Wilcoxon.test p value=",format(w2$p.value, scientific = TRUE)),ylab="log2 TMM-normalized CPM ",
            names=c("IDC(n=7)","ILC(n=14)"),col=c("darkblue","darkred"),notch=FALSE)
  })
  
  output$gi <-renderPlot({
    gene<-input$gene
    chosedGene <- subset(PanMets, (external_gene_name %in% gene))
    
    mdat = melt(chosedGene, id.vars=("external_gene_name"))
    
    mdat_wkey<- merge(mdat,Key_file,by.x = "variable",by.y = "ID")
    colnames(mdat_wkey)[3] <- "log2CPM"
    
    mdat_tumors <- mdat_wkey
    #PAIRS
    mdat_pairs <- mdat_tumors[!(is.na(mdat_tumors$pair) | mdat_tumors$pair==""),]
    mdat_pairs <- mdat_pairs[!(mdat_pairs$subtype=="unknown" |  mdat_pairs$pair==""),]
    mdat_pairs$group <- factor(mdat_pairs$group, c("Primary", "MET")) #change order so primary appear left on plot
    mdat_pairs$site <- factor(mdat_pairs$site, c("Brain","Bone", "Ovary","Breast","GI"))
    mdat_pairs$ER_status <- factor(mdat_pairs$ER_status, c("Pos","Neg","Unk"))
    GI <- mdat_pairs[mdat_pairs$code== "Gi",]
    test<- subset(GI, (external_gene_name %in% gene))
    test.primary <- subset(test, (group %in% "Primary"))
    test.met <- subset(test, (group %in% "MET"))
    w<-wilcox.test(test.primary$log2CPM, test.met$log2CPM,paired = TRUE, alternative = "two.sided")
    
    ggplot(subset(GI, (external_gene_name %in% gene)), aes(x=group, y=log2CPM, color= subtype,group=pair))  +
      geom_point(size=3) +geom_line(aes(linetype=ER_status),size=0.5) + theme(axis.text=element_text(size=15,face="bold"),
                                                                              axis.title=element_text(size=14,face="bold")) +
      ggtitle(paste(gene,"Expression in Paired GI Metastases")) + 
      theme(plot.title = element_text(color="black", face="bold", size=14, hjust=.5,vjust=1)) +
      labs(x="",y="log2(TMM-normalized CPM +1)",caption = paste("Paired Wilcoxon.test p value=",format(w$p.value, scientific = TRUE)))+
      scale_color_manual(values=c("purple4", "red4")) +
      scale_linetype_manual(values=c(1,3))+
      theme(legend.title=element_text(size=14,face="bold"))+
      theme(legend.text=element_text(size=12,face="bold"))
  })
  
  output$br <-renderPlot({
    gene<-input$gene
    chosedGene <- subset(PanMets, (external_gene_name %in% gene))
    
    mdat = melt(chosedGene, id.vars=("external_gene_name"))
    
    mdat_wkey<- merge(mdat,Key_file,by.x = "variable",by.y = "ID")
    colnames(mdat_wkey)[3] <- "log2CPM"
    
    mdat_tumors <- mdat_wkey
    #PAIRS
    mdat_pairs <- mdat_tumors[!(is.na(mdat_tumors$pair) | mdat_tumors$pair==""),]
    mdat_pairs <- mdat_pairs[!(mdat_pairs$subtype=="unknown" |  mdat_pairs$pair==""),]
    mdat_pairs$group <- factor(mdat_pairs$group, c("Primary", "MET")) #change order so primary appear left on plot
    mdat_pairs$site <- factor(mdat_pairs$site, c("Brain","Bone", "Ovary","Breast","GI"))
    mdat_pairs$ER_status <- factor(mdat_pairs$ER_status, c("Pos","Neg","Unk"))
    Brain <- mdat_pairs[mdat_pairs$code== "Br",]
    test<- subset(Brain, (external_gene_name %in% gene))
    test.primary <- subset(test, (group %in% "Primary"))
    test.met <- subset(test, (group %in% "MET"))
    w<-wilcox.test(test.primary$log2CPM, test.met$log2CPM,paired = TRUE, alternative = "two.sided")
    
    ggplot(subset(Brain, (external_gene_name %in% gene)), aes(x=group, y=log2CPM, color= subtype,group=pair))  +
      geom_point(size=3) +geom_line(aes(linetype=ER_status),size=0.5) +theme(axis.text=element_text(size=15,face="bold"),
                                                                             axis.title=element_text(size=14,face="bold")) +
      ggtitle(paste(gene,"Expression in Paired Brain Metastases")) + 
      theme(plot.title = element_text(color="black", face="bold", size=14, hjust=.5,vjust=1)) +
      labs(x="",y="log2(TMM-normalized CPM +1)",caption = paste("Paired Wilcoxon.test p value=",format(w$p.value, scientific = TRUE)))+
      scale_color_manual(values=c("purple4", "red4")) +
      scale_linetype_manual(values=c(1,3)) +
      theme(legend.title=element_text(size=14,face="bold"))+
      theme(legend.text=element_text(size=12,face="bold"))
  })
  
  output$br.er <-renderPlot({
    gene<-input$gene
    chosedGene <- subset(PanMets, (external_gene_name %in% gene))
    
    mdat = melt(chosedGene, id.vars=("external_gene_name"))
    
    mdat_wkey<- merge(mdat,Key_file,by.x = "variable",by.y = "ID")
    colnames(mdat_wkey)[3] <- "log2CPM"
    
    mdat_tumors <- mdat_wkey
    #PAIRS
    mdat_pairs <- mdat_tumors[!(is.na(mdat_tumors$pair) | mdat_tumors$pair==""),]
    mdat_pairs <- mdat_pairs[!(mdat_pairs$subtype=="unknown" |  mdat_pairs$pair==""),]
    mdat_pairs$group <- factor(mdat_pairs$group, c("Primary", "MET")) #change order so primary appear left on plot
    mdat_pairs$site <- factor(mdat_pairs$site, c("Brain","Bone", "Ovary","Breast","GI"))
    mdat_pairs$ER_status <- factor(mdat_pairs$ER_status, c("Pos","Neg","Unk"))
    Brain <- mdat_pairs[mdat_pairs$code== "Br",]
    Brain.MET <- Brain[Brain$group== "MET",]
    Brain.MET.gene <- subset(Brain.MET, external_gene_name %in% gene)
    Brain.MET.ERpos <- Brain.MET.gene[Brain.MET.gene$ER_status=="Pos",]
    Brain.MET.ERneg <- Brain.MET.gene[Brain.MET.gene$ER_status=="Neg",]
    w3<-wilcox.test(Brain.MET.ERpos$log2CPM,Brain.MET.ERneg$log2CPM,paired = FALSE, alternative = "two.sided")
    
    boxplot(Brain.MET.ERpos$log2CPM,Brain.MET.ERneg$log2CPM,main= paste(gene,"expression in Brain Mets(ERpos vs ERneg)"),cex.main=1.5,outline = TRUE,pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5),
            log = "",
            xlab=paste("Non-Paired Wilcoxon.test p value=",format(w3$p.value, scientific = TRUE)),ylab="log2 TMM-normalized CPM ",
            names=c("ERpos(n=9)","ERneg(n=13)"),col=c("brown3","cyan4"),notch=FALSE)
  })
  
  output$bo <-renderPlot({
    gene<-input$gene
    chosedGene <- subset(PanMets, (external_gene_name %in% gene))
    
    mdat = melt(chosedGene, id.vars=("external_gene_name"))
    
    mdat_wkey<- merge(mdat,Key_file,by.x = "variable",by.y = "ID")
    colnames(mdat_wkey)[3] <- "log2CPM"
    
    mdat_tumors <- mdat_wkey
    #PAIRS
    mdat_pairs <- mdat_tumors[!(is.na(mdat_tumors$pair) | mdat_tumors$pair==""),]
    mdat_pairs <- mdat_pairs[!(mdat_pairs$subtype=="unknown" |  mdat_pairs$pair==""),]
    mdat_pairs$group <- factor(mdat_pairs$group, c("Primary", "MET")) #change order so primary appear left on plot
    mdat_pairs$site <- factor(mdat_pairs$site, c("Brain","Bone", "Ovary","Breast","GI"))
    mdat_pairs$ER_status <- factor(mdat_pairs$ER_status, c("Pos","Neg","Unk"))
    Bone <- mdat_pairs[mdat_pairs$code== "Bo",]
    test<- subset(Bone, (external_gene_name %in% gene))
    test.primary <- subset(test, (group %in% "Primary"))
    test.met <- subset(test, (group %in% "MET"))
    w<-wilcox.test(test.primary$log2CPM, test.met$log2CPM,paired = TRUE, alternative = "two.sided")
    
    ggplot(subset(Bone, (external_gene_name %in% gene)), aes(x=group, y=log2CPM, color= subtype,group=pair))  +
      geom_point(size=3) +geom_line(aes(linetype=ER_status),size=0.5) + theme(axis.text=element_text(size=15,face="bold"),
                                                                              axis.title=element_text(size=14,face="bold")) +
      ggtitle(paste(gene,"Expression in Paired Bone Metastases")) + 
      theme(plot.title = element_text(color="black", face="bold", size=14, hjust=.5,vjust=1)) +
      labs(x="",y="log2(TMM-normalized CPM +1)",caption = paste("Paired Wilcoxon.test p value=",format(w$p.value, scientific = TRUE))) +
      scale_color_manual(values=c("purple4", "red4","darkorange")) +
      scale_linetype_manual(values=c(1,3))+
      theme(legend.title=element_text(size=14,face="bold"))+
      theme(legend.text=element_text(size=12,face="bold"))
  })
})