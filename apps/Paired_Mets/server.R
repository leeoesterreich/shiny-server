library(shiny)
library(ggplot2)
library(data.table)
library(reshape)

load("PanMets_BrBo.Rda")
load("PanMets_OvGi.Rda")
Key_file= read.csv("SampleAnnotation.csv")

shinyServer(function(input, output) {
  
  #all
  output$all <-renderPlot({
    gene<-input$gene
    chosedGene_og <- subset(PanMets_OvGi, (external_gene_name %in% gene))
    mdat_og = melt(chosedGene_og, id.vars=("external_gene_name"))
    chosedGene_bb <- subset(PanMets_BrBo, (external_gene_name %in% gene))
    mdat_bb = melt(chosedGene_bb, id.vars=("external_gene_name"))
    mdat<-rbind(mdat_og,mdat_bb)
    
    mdat_wkey<- merge(mdat,Key_file,by.x = "variable",by.y = "ID")
    colnames(mdat_wkey)[3] <- "log2CPM"
    
    mdat_tumors <- mdat_wkey
    #PAIRS
    mdat_pairs <- mdat_tumors[!(is.na(mdat_tumors$pair) | mdat_tumors$pair==""),]
    mdat_pairs$group <- factor(mdat_pairs$group, c("Primary", "MET")) #change order so primary appear left on plot
    mdat_pairs$site <- factor(mdat_pairs$site, c("Brain","Bone", "Ovary","Breast","GI"))
    mdat_pairs$metSite <- ifelse(mdat_pairs$code=="Bo", "Bone", ifelse(mdat_pairs$code=="Br", "Brain",ifelse(mdat_pairs$code=="Gi", "GI", ifelse(mdat_pairs$code=="Ov", "Ovary", NA))))
    mdat_pairs$ER_status <- factor(mdat_pairs$ER_status, c("Pos","Neg","Unk"))
    
    if(input$erInput=="ER+") mdat_pairs<-subset(mdat_pairs, ER_status=="Pos") 
    if(input$erInput=="ER-") mdat_pairs<-subset(mdat_pairs, ER_status=="Neg") 
    
    test<- subset(mdat_pairs, (external_gene_name %in% gene))
    test.primary <- subset(test, (group %in% "Primary"))
    test.met <- subset(test, (group %in% "MET"))
    w<-wilcox.test(test.primary$log2CPM, test.met$log2CPM,paired = TRUE, alternative = "two.sided")
    
    ggplot(subset(mdat_pairs, (external_gene_name %in% gene)), aes(x=group, y=log2CPM, color=metSite,group=pair))  +
      geom_point(size=3) +geom_line(aes(linetype=ER_status),size=0.5) + theme(axis.text=element_text(size=15,face="bold"),
                                                                              axis.title=element_text(size=14,face="bold")) +
      ggtitle(paste(gene,"Expression in Paired Metastases")) + 
      theme(plot.title = element_text(color="black", face="bold", size=14, hjust=.5,vjust=1)) +
      labs(x=paste("paired Wilcoxon.test p-value=",format(w$p.value, digits = 2)),y="log2(TMM-normalized CPM +1)") +
      scale_linetype_manual(values= c(1,3)) +
      theme(legend.title=element_text(size=14,face="bold"))+
      theme(legend.text=element_text(size=12,face="bold"))+
      theme(axis.title=element_text(size=12,face="plain"))
  })
  
  #Brain
  output$br <-renderPlot({
    gene<-input$gene
    chosedGene <- subset(PanMets_BrBo, (external_gene_name %in% gene))
    
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
    
    if(input$erInput=="ER+") mdat_pairs<-subset(mdat_pairs, ER_status=="Pos") 
    if(input$erInput=="ER-") mdat_pairs<-subset(mdat_pairs, ER_status=="Neg") 
    
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
      labs(x=paste("paired Wilcoxon.test p-value=",format(w$p.value, digits = 2)),y="log2(TMM-normalized CPM +1)") +
      scale_color_manual(values=c("purple4", "red4")) +
      scale_linetype_manual(values=c(1,3)) +
      theme(legend.title=element_text(size=14,face="bold"))+
      theme(legend.text=element_text(size=12,face="bold"))+
      theme(axis.title=element_text(size=12,face="plain"))
  })
  
  #Bone
  output$bo <-renderPlot({
    gene<-input$gene
    chosedGene <- subset(PanMets_BrBo, (external_gene_name %in% gene))
    
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
    
    if(input$erInput=="ER+") mdat_pairs<-subset(mdat_pairs, ER_status=="Pos") 
    
    #No ER- bone mets
    if(input$erInput!="ER-"){
      
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
        labs(x=paste("paired Wilcoxon.test p-value=",format(w$p.value, digits = 2)),y="log2(TMM-normalized CPM +1)") +
        scale_color_manual(values=c("purple4", "red4","darkorange")) +
        scale_linetype_manual(values=c(1,3))+
        theme(legend.title=element_text(size=14,face="bold"))+
        theme(legend.text=element_text(size=12,face="bold"))+
        theme(axis.title=element_text(size=12,face="plain"))
    }
    
    #if(input$erInput=="ER-"){
    #  
    #}
  })
  
  #Ovary
  output$ov <-renderPlot({
    gene<-input$gene
    chosedGene <- subset(PanMets_OvGi, (external_gene_name %in% gene))
    
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
    
    if(input$erInput=="ER+") mdat_pairs<-subset(mdat_pairs, ER_status=="Pos") 
    if(input$erInput=="ER-") mdat_pairs<-subset(mdat_pairs, ER_status=="Neg") 
    
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
      labs(x=paste("paired Wilcoxon.test p-value=",format(w$p.value, digits = 2)),y="log2(TMM-normalized CPM +1)") +
      scale_color_manual(values= c("purple4", "red4","darkorange")) +
      scale_linetype_manual(values= c(1,3)) +
      theme(legend.title=element_text(size=14,face="bold"))+
      theme(legend.text=element_text(size=12,face="bold"))+
      theme(axis.title=element_text(size=12,face="plain"))
  })
  
  #GI
  output$gi <-renderPlot({
    gene<-input$gene
    chosedGene <- subset(PanMets_OvGi, (external_gene_name %in% gene))
    
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
    
    if(input$erInput=="ER+") mdat_pairs<-subset(mdat_pairs, ER_status=="Pos") 
    if(input$erInput=="ER-") mdat_pairs<-subset(mdat_pairs, ER_status=="Neg") 
    
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
      labs(x=paste("paired Wilcoxon.test p-value=",format(w$p.value, digits = 2)),y="log2(TMM-normalized CPM +1)") +
      scale_color_manual(values=c("purple4", "red4")) +
      scale_linetype_manual(values=c(1,3))+
      theme(legend.title=element_text(size=14,face="bold"))+
      theme(legend.text=element_text(size=12,face="bold"))+
      theme(axis.title=element_text(size=12,face="plain"))
  })
  
  output$br.er <-renderPlot({
    gene<-input$gene
    chosedGene <- subset(PanMets_BrBo, (external_gene_name %in% gene))
    
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
    
    boxplot(Brain.MET.ERpos$log2CPM,Brain.MET.ERneg$log2CPM,main= paste(gene,"expression in\n Brain Mets"),cex.main=1,outline = TRUE,pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5),
            log = "",
            xlab=paste("Non-Paired Wilcoxon.test\np-value=",format(w3$p.value, digits=2)),ylab="log2 TMM-normalized CPM ",
            names=c("ERpos\n(n=9)","ERneg\n(n=13)"),col=c("brown3","cyan4"),notch=FALSE)
  })
  
  output$ov.sb <-renderPlot({
    gene<-input$gene
    chosedGene <- subset(PanMets_OvGi, (external_gene_name %in% gene))
    
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
    
    boxplot(ovaries.MET.IDC$log2CPM,ovaries.MET.ILC$log2CPM,main= paste(gene,"expression in\n Ovarian Mets"),cex.main=1,outline = TRUE,pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5),
            log = "",
            xlab=paste("Non-Paired Wilcoxon.test\np-value=",format(w2$p.value, digits = 2)),ylab="log2 TMM-normalized CPM ",
            names=c("IDC\n(n=7)","ILC\n(n=14)"),col=c("darkblue","darkred"),notch=FALSE)
  })
  
})