load("primary_tpm.Rda")
load("normal_tpm.Rda")
load("bind.Rda")
#group<-as.list(colnames(df[,23369:23372]))
shinyServer(function(input, output) {
  
  PAM<-reactive({
    BC<-df[,c(input$gene,"pam50")]
    BC$pam50<-as.factor(as.character(BC$pam50))
    return(BC)
  })
  
  DataPrimaries_menu1<-reactive({
    BC<-df[,c(input$gene,"er_status_by_ihc","Histology_textCode","pam50")]
    BC<-cbind(rownames(BC),BC)
    colnames(BC)<-c("Patient_ID (primaries)",input$gene,"er_status_by_ihc","Histology_textCode","pam50")
    if(input$dataset=="Breast Cancer"){
      return(BC)
    }else if(input$dataset=="ER+"){
      sele<-BC[BC$er_status_by_ihc=="Positive",]
      sele$er_status_by_ihc<-as.factor(as.character(sele$er_status_by_ihc))
      return(sele)
    }else if(input$dataset=="ER-"){
      sele<-BC[BC$er_status_by_ihc=="Negative",]
      sele$er_status_by_ihc<-as.factor(as.character(sele$er_status_by_ihc))
      return(sele)
    }else if(input$dataset=="ILC"){
      sele<-BC[BC$Histology_textCode=="ILC",]
      sele<-na.omit(sele)
      sele$Histology_textCode<-as.factor(as.character(sele$Histology_textCode))
      return(sele)
    }else if(input$dataset=="IDC"){
      sele<-BC[BC$Histology_textCode=="IDC",]
      sele<-na.omit(sele)
      sele$Histology_textCode<-as.factor(as.character(sele$Histology_textCode))
      return(sele)
    }else if(input$dataset=="LumA"){
      sele<-BC[BC$pam50=="LumA",]
      sele$pam50<-as.factor(as.character(sele$pam50))
      return(sele)
    }else if(input$dataset=="LumB"){
      sele<-BC[BC$pam50=="LumB",]
      sele$pam50<-as.factor(as.character(sele$pam50))
      return(sele)
    }else if(input$dataset=="Her2"){
      sele<-BC[BC$pam50=="Her2",]
      sele$pam50<-as.factor(as.character(sele$pam50))
      return(sele)
    }else if(input$dataset=="Basal"){
      sele<-BC[BC$pam50=="Basal",]
      sele$pam50<-as.factor(as.character(sele$pam50))
      return(sele)
    }
  })
  
  Primaries<-reactive({  # for statistical test
    data<-DataPrimaries_menu1()
    if(input$group=="ER+ vs ER-"){
      sele<-data[data$er_status_by_ihc=="Positive"|data$er_status_by_ihc=="Negative",]
      sele$er_status_by_ihc<-as.factor(as.character(sele$er_status_by_ihc))
      return(sele)
    }else if(input$group=="ILC vs IDC"){
      sele<-data[data$Histology_textCode=="ILC"|data$Histology_textCode=="IDC",]
      sele$Histology_textCode<-as.factor(as.character(sele$Histology_textCode))
      sele<-na.omit(sele)
      return(sele)
    }else if(input$group=="LumA vs LumB"){
      sele<-data[data$pam50=="LumA"|data$pam50=="LumB",]
      sele$pam50<-as.factor(as.character(sele$pam50))
      return(sele)
    }
  })
  
  NT<-reactive({ # tumor vs. adj.tissue
    BC<-bind[,c(input$gene,"er_status_by_ihc")]
    BC<-cbind(rownames(BC),BC)
    colnames(BC)<-c("Patient_ID",input$gene,"er_status_by_ihc")
    sele<-BC[BC$er_status_by_ihc=="Tumor Adjacent Normal Tissue"|BC$er_status_by_ihc=="T",]
    return(sele)
  })
  
  DataBind_menu1<-reactive({
    BC<-bind[,c(input$gene,"er_status_by_ihc","Histology_textCode","pam50")]
    BC<-cbind(rownames(BC),BC)
    colnames(BC)<-c("Patient_ID (primaries)",input$gene,"er_status_by_ihc","Histology_textCode","pam50")
    
    if(input$dataset=="Breast Cancer"){
      return(BC)
    }else if(input$dataset=="ER+"){
      sele<-BC[BC$er_status_by_ihc=="Positive"|BC$er_status_by_ihc=="Tumor Adjacent Normal Tissue"|BC$er_status_by_ihc=="T",]
      sele$er_status_by_ihc<-as.factor(as.character(sele$er_status_by_ihc))
      return(sele)
    }else if(input$dataset=="ER-"){
      sele<-BC[BC$er_status_by_ihc=="Negative"|BC$er_status_by_ihc=="Tumor Adjacent Normal Tissue"|BC$er_status_by_ihc=="T",]
      sele$er_status_by_ihc<-as.factor(as.character(sele$er_status_by_ihc))
      return(sele)
    }else if(input$dataset=="ILC"){
      sele<-BC[BC$Histology_textCode=="ILC"|BC$Histology_textCode=="Tumor Adjacent Normal Tissue"|BC$Histology_textCode=="T",]
      sele<-na.omit(sele)
      sele$Histology_textCode<-as.factor(as.character(sele$Histology_textCode))
      return(sele)
    }else if(input$dataset=="IDC"){
      sele<-BC[BC$Histology_textCode=="IDC"|BC$Histology_textCode=="Tumor Adjacent Normal Tissue"|BC$Histology_textCode=="T",]
      sele<-na.omit(sele)
      sele$Histology_textCode<-as.factor(as.character(sele$Histology_textCode))
      return(sele)
    }else if(input$dataset=="LumA"){
      sele<-BC[BC$pam50=="LumA"|BC$pam50=="Tumor Adjacent Normal Tissue"|BC$pam50=="T",]
      sele$pam50<-as.factor(as.character(sele$pam50))
      return(sele)
    }else if(input$dataset=="LumB"){
      sele<-BC[BC$pam50=="LumB"|BC$pam50=="Tumor Adjacent Normal Tissue"|BC$pam50=="T",]
      sele$pam50<-as.factor(as.character(sele$pam50))
      return(sele)
    }else if(input$dataset=="Her2"){
      sele<-BC[BC$pam50=="Her2"|BC$pam50=="Tumor Adjacent Normal Tissue"|BC$pam50=="T",]
      sele$pam50<-as.factor(as.character(sele$pam50))
      return(sele)
    }else if(input$dataset=="Basal"){
      sele<-BC[BC$pam50=="Basal"|BC$pam50=="Tumor Adjacent Normal Tissue"|BC$pam50=="T",]
      sele$pam50<-as.factor(as.character(sele$pam50))
      return(sele)
    }
  })
  
  Bind<-reactive({ # for boxplot
    data<-DataBind_menu1()
    if(input$group=="ER+ vs ER-"){
      sele<-data[data$er_status_by_ihc=="Positive"|data$er_status_by_ihc=="Negative"|data$er_status_by_ihc=="Tumor Adjacent Normal Tissue"|data$er_status_by_ihc=="T",]
      sele$er_status_by_ihc<-as.factor(as.character(sele$er_status_by_ihc))
      return(sele)
    }else if(input$group=="ILC vs IDC"){
      sele<-data[data$Histology_textCode=="ILC"|data$Histology_textCode=="IDC"|data$Histology_textCode=="Tumor Adjacent Normal Tissue"|data$Histology_textCode=="T",]
      sele$Histology_textCode<-as.factor(as.character(sele$Histology_textCode))
      sele<-na.omit(sele)
      return(sele)
    }else if(input$group=="LumA vs LumB"){
      sele<-data[data$pam50=="LumA"|data$pam50=="LumB"|data$pam50=="Tumor Adjacent Normal Tissue"|data$pam50=="T",]
      sele$pam50<-as.factor(as.character(sele$pam50))
      return(sele)
    }
  })
  
  DataNormal_menu1<-reactive({
    BC<-normal_tpm[,c(input$gene,"er_status_by_ihc","Histology_textCode","pam50")]
    BC<-cbind(rownames(BC),BC)
    colnames(BC)<-c("Patient_ID",input$gene,"er_status_by_ihc","Histology_textCode","pam50")
    if(input$dataset=="Breast Cancer"){
      return(BC)
    }else if(input$dataset=="ER+"){
      sele<-BC[BC$er_status_by_ihc=="Positive",]
      sele$er_status_by_ihc<-as.factor(as.character(sele$er_status_by_ihc))
      return(sele)
    }else if(input$dataset=="ER-"){
      sele<-BC[BC$er_status_by_ihc=="Negative",]
      sele$er_status_by_ihc<-as.factor(as.character(sele$er_status_by_ihc))
      return(sele)
    }else if(input$dataset=="ILC"){
      sele<-BC[BC$Histology_textCode=="ILC",]
      sele<-na.omit(sele)
      sele$Histology_textCode<-as.factor(as.character(sele$Histology_textCode))
      return(sele)
    }else if(input$dataset=="IDC"){
      sele<-BC[BC$Histology_textCode=="IDC",]
      sele<-na.omit(sele)
      sele$Histology_textCode<-as.factor(as.character(sele$Histology_textCode))
      return(sele)
    }else if(input$dataset=="LumA"){
      sele<-BC[BC$pam50=="LumA",]
      sele$pam50<-as.factor(as.character(sele$pam50))
      return(sele)
    }else if(input$dataset=="LumB"){
      sele<-BC[BC$pam50=="LumB",]
      sele$pam50<-as.factor(as.character(sele$pam50))
      return(sele)
    }else if(input$dataset=="Her2"){
      sele<-BC[BC$pam50=="Her2",]
      sele$pam50<-as.factor(as.character(sele$pam50))
      return(sele)
    }else if(input$dataset=="Basal"){
      sele<-BC[BC$pam50=="Basal",]
      sele$pam50<-as.factor(as.character(sele$pam50))
      return(sele)
    }
  })
  
  Normal<-reactive({
    data<-DataNormal_menu1()
    if(input$group=="ER+ vs ER-"){
      sele<-data[data$er_status_by_ihc=="Positive"|data$er_status_by_ihc=="Negative",]
      sele$er_status_by_ihc<-as.factor(as.character(sele$er_status_by_ihc))
      return(sele)
    }else if(input$group=="ILC vs IDC"){
      sele<-data[data$Histology_textCode=="ILC"|data$Histology_textCode=="IDC",]
      sele$Histology_textCode<-as.factor(as.character(sele$Histology_textCode))
      sele<-na.omit(sele)
      return(sele)
    }else if(input$group=="LumA vs LumB"){
      sele<-data[data$pam50=="LumA"|data$pam50=="LumB",]
      sele$pam50<-as.factor(as.character(sele$pam50))
      return(sele)
    }
  })
 
    output$view1<-renderDataTable({
      Primaries()
    }, options = list(lengthMenu = c(20,50,100), pageLength=20)
    )
    
    output$plot1<-renderPlot({
      form<-Bind()
      if(input$group=="ER+ vs ER-"){
        a<-length(form[form$er_status_by_ihc=="Negative",3])
        b<-length(form[form$er_status_by_ihc=="Positive",3])
        par(mgp = c(3, 1.5, 0))
        boxplot(Bind()[,2]~Bind()[,3],
                main=paste(input$gene,"Expression in Primary Tumor",sep=" "),
                ylab="log2 (TPM+1)",
                col=ifelse(levels(Bind()$er_status_by_ihc)=="Tumor Adjacent Normal Tissue","grey90",rgb(0.8,0.1,0.3,0.6)),
                names=c(paste("ER-\n(n=",a,")",sep=""),paste("ER+\n(n=",b,")",sep=""),"all Tumors\n(n=1095)","all Adj. Normal\n(n=113)")
                )
      }else if(input$group=="ILC vs IDC"){
        a<-length(form[form$Histology_textCode=="IDC",4])
        b<-length(form[form$Histology_textCode=="ILC",4])
        par(mgp = c(3, 1.5, 0))
        boxplot(Bind()[,2]~Bind()[,4],
                main=paste(input$gene,"Expression in Primary Tumor",sep=" "),  
                ylab="log2 (TPM+1)",
                col=ifelse(levels(Bind()$Histology_textCode)=="Tumor Adjacent Normal Tissue","grey90",rgb(0.8,0.1,0.3,0.6)),
                names=c(paste("IDC\n(n=",a,")",sep=""), paste("ILC\n(n=",b,")",sep=""),"all Tumors\n(n=1095)","all Adj. Normal\n(n=113)")
                )
      }else if(input$group=="LumA vs LumB"){
        a<-length(form[form$pam50=="LumA",5])
        b<-length(form[form$pam50=="LumB",5])
        par(mgp = c(3, 1.5, 0))
        boxplot(Bind()[,2]~Bind()[,5],
                main=paste(input$gene,"Expression in Primary Tumor",sep=" "),  
                ylab="log2 (TPM+1)",
                col=ifelse(levels(Bind()$pam50)=="Tumor Adjacent Normal Tissue","grey90",rgb(0.8,0.1,0.3,0.6)),
                names=c(paste("LumA\n(n=",a,")",sep=""), paste("LumB\n(n=",b,")",sep=""),"all Tumors\n(n=1095)","all Adj. Normal\n(n=113)")
                )
      }
    })
    
    output$plot2<-renderPlot({
      form<-Normal()
      if(input$group=="ER+ vs ER-"){
        a<-length(form[form$er_status_by_ihc=="Negative",3])
        b<-length(form[form$er_status_by_ihc=="Positive",3])
        par(mgp = c(3, 1.5, 0))
        boxplot(Normal()[,2]~Normal()[,3],
                main=paste(input$gene,"Expression in Adjacent Tissue",sep=" "),
                ylab="log2 (TPM+1)",
                names=c(paste("ER-\n(n=",a,")",sep=""), paste("ER+\n(n=",b,")",sep="")),
                col="grey90"
                )
      }else if(input$group=="ILC vs IDC"){
        a<-length(form[form$Histology_textCode=="IDC",4])
        b<-length(form[form$Histology_textCode=="ILC",4])
        par(mgp = c(3, 1.5, 0))
        boxplot(Normal()[,2]~Normal()[,4],
                main=paste(input$gene,"Expression in Adjacent Tissue",sep=" "),  
                ylab="log2 (TPM+1)",
                names=c(paste("IDC\n(n=",a,")",sep=""), paste("ILC\n(n=",b,")",sep="")),
                col="grey90"
                )
      }else if(input$group=="LumA vs LumB"){
        a<-length(form[form$pam50=="LumA",5])
        b<-length(form[form$pam50=="LumB",5])
        par(mgp = c(3, 1.5, 0))
        boxplot(Normal()[,2]~Normal()[,5],
                main=paste(input$gene,"Expression in Adjacent Tissue",sep=" "),  
                ylab="log2 (TPM+1)",
                names=c(paste("LumA\n(n=",a,")",sep=""), paste("LumB\n(n=",b,")",sep="")),
                col="grey90"
                )
      }
    })
    
    output$plot3<-renderPlot({
      form<-PAM()
      a<-length(form[form$pam50=="Basal",2])
      b<-length(form[form$pam50=="Her2",2])
      c<-length(form[form$pam50=="LumA",2])
      d<-length(form[form$pam50=="LumB",2])
      e<-length(form[form$pam50=="Normal",2])
      par(mgp = c(3, 1.5, 0))
      boxplot(PAM()[,1]~PAM()[,2],
              main=paste(input$gene,"Expression in all Primary Tumors (PAM50)",sep=" "),
              ylab="log2 (TPM+1)",
              names=c(paste("Basal\n(n=",a,")",sep=""),paste("Her2\n(n=",b,")",sep=""),paste("LumA\n(n=",c,")",sep=""), 
                      paste("LumB\n(n=",d,")",sep=""),paste("Normal-like\n(n=",e,")",sep="")),
              col=rgb(0.8,0.1,0.3,0.6)
              )
    })
    
    output$state<-renderPrint({
      cat("TCGA data downloaded from:\nhttps://www.ncbi.nlm.nih.gov/pubmed/26209429") 
    })
    
    output$stats1<-renderPrint({
      if(input$group=="ER+ vs ER-"){
        cat(paste(input$gene,"Expression in Primary Tumor:\nER+ vs ER- in",input$dataset,"patients",sep=" "),
            paste("\nMann-Whitney U test (ER+ vs ER-): p-value=",signif(wilcox.test(Primaries()[,2]~Primaries()[,3])$p.value,digits=3),sep=""),
            paste("\nSummary of ER-",input$dataset,"[log2(TPM+1)]:",sep=" "),
            paste("Mean =",round(mean(Primaries()[Primaries()$er_status_by_ihc=="Negative",2]),digits=3),sep=" "),
            paste("Median =",round(median(Primaries()[Primaries()$er_status_by_ihc=="Negative",2]),digits=3),sep=" "),
            paste("\nSummary of ER+",input$dataset,"[log2(TPM+1)]:",sep=" "),
            paste("Mean =",round(mean(Primaries()[Primaries()$er_status_by_ihc=="Positive",2]),digits=3),sep=" "),
            paste("Median =",round(median(Primaries()[Primaries()$er_status_by_ihc=="Positive",2]),digits=3),sep=" "),
            sep = "\n")
      }else if(input$group=="ILC vs IDC"){
        cat(paste(input$gene,"Expression in Primary Tumor:\nILC vs IDC in",input$dataset,"patients",sep=" "),
            paste("\nMann-Whitney U test (ILC vs IDC): p-value=",signif(wilcox.test(Primaries()[,2]~Primaries()[,4])$p.value,digits=3),sep=""),
            paste("\nSummary of IDC",input$dataset,"[log2(TPM+1)]:",sep=" "),
            paste("Mean =",round(mean(Primaries()[Primaries()$Histology_textCode=="IDC",2]),digits=3),sep=" "),
            paste("Median =",round(median(Primaries()[Primaries()$Histology_textCode=="IDC",2]),digits=3),sep=" "),
            paste("\nSummary of ILC",input$dataset,"[log2(TPM+1)]:",sep=" "),
            paste("Mean =",round(mean(Primaries()[Primaries()$Histology_textCode=="ILC",2]),digits=3),sep=" "),
            paste("Median =",round(median(Primaries()[Primaries()$Histology_textCode=="ILC",2]),digits=3),sep=" "),
            sep = "\n")
      }else if(input$group=="LumA vs LumB"){
        cat(paste(input$gene,"Expression in Primary Tumor:\nLumA vs LumB in",input$dataset,"patients",sep=" "),
            paste("\nMann-Whitney U test (LumA vs LumB): p-value=",signif(wilcox.test(Primaries()[,2]~Primaries()[,5])$p.value,digits=3),sep=""),
            paste("\nSummary of LumA",input$dataset,"[log2(TPM+1)]:",sep=" "),
            paste("Mean =",round(mean(Primaries()[Primaries()$pam50=="LumA",2]),digits=3),sep=" "),
            paste("Median =",round(median(Primaries()[Primaries()$pam50=="LumA",2]),digits=3),sep=" "),
            paste("\nSummary of LumB",input$dataset,"[log2(TPM+1)]:",sep=" "),
            paste("Mean =",round(mean(Primaries()[Primaries()$pam50=="LumB",2]),digits=3),sep=" "),
            paste("Median =",round(median(Primaries()[Primaries()$pam50=="LumB",2]),digits=3),sep=" "),
            sep = "\n")
      }
    })
    
    output$stats2<-renderPrint({
      if(input$group=="ER+ vs ER-"){
        cat(paste(input$gene,"Expression in Adjacent Tissue:\nER+ vs ER- in",input$dataset,"patients",sep=" "),
            paste("\nMann-Whitney U test (ER+ vs ER-): p-value=",signif(wilcox.test(Normal()[,2]~Normal()[,3])$p.value,digits=3),sep=""),
            sep="\n")
      }else if(input$group=="ILC vs IDC"){
        cat(paste(input$gene,"Expression in Adjacent Tissue:\nILC vs IDC in",input$dataset,"patients",sep=" "),
            paste("\nMann-Whitney U test (ILC vs IDC): p-value=",signif(wilcox.test(Normal()[,2]~Normal()[,4])$p.value,digits=3),sep=""),
            sep="\n")
      }else if(input$group=="LumA vs LumB"){
        cat(paste(input$gene,"Expression in Adjacent Tissue:\nLumA vs LumB in",input$dataset,"patients",sep=" "),
            paste("\nMann-Whitney U test (LumA vs LumB): p-value=",signif(wilcox.test(Normal()[,2]~Normal()[,5])$p.value,digits=3),sep=""),
            sep="\n")
      }
    })
    
    output$stats3<-renderPrint({
      cat(paste(input$gene,"Expression in Tumor vs Adjacent Tissue"),
          paste("\nMann-Whitney U test: p-value=",signif(wilcox.test(NT()[,2]~NT()[,3])$p.value,digits=3),sep=""),
          "\nSummary of Tumor [log2(TPM+1)]:",
          paste("Mean =",round(mean(NT()[NT()$er_status_by_ihc=="T",2]),digits=3),sep=" "),
          paste("Median =",round(median(NT()[NT()$er_status_by_ihc=="T",2]),digits=3),sep=" "),
          "\nSummary of Tumor Adjacent Tissue [log2(TPM+1)]:",
          paste("Mean =",round(mean(NT()[NT()$er_status_by_ihc=="Tumor Adjacent Normal Tissue",2]),digits=3),sep=" "),
          paste("Median =",round(median(NT()[NT()$er_status_by_ihc=="Tumor Adjacent Normal Tissue",2]),digits=3),sep=" "),
          sep="\n")
    })
    
    output$stats4<-renderPrint({
      cat(paste(input$gene,"Expression in all Primary Tumors",sep=" "),
          "Tukey multiple comparisons of means",
          "95% family-wise confidence level (p-value)\n",
          sep="\n")
      a1<-aov(PAM()[,1]~PAM()[,2])
      tt<-TukeyHSD(a1)
      result<-data.frame(tt$`PAM()[, 2]`)
      result["p.adj"]
    })
    
  }
)