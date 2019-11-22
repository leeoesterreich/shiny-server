library(shiny)

shinyUI(fluidPage(
  titlePanel("Gene expression from RNA-seq of different ESR1 mutation cell models"),
  
  sidebarLayout(
    sidebarPanel(
      
      textInput("gene", label= "Enter Gene Symbol in CAPITAL letters", value = "ESR1"),

      helpText("RNA-seq data were extracted from three resources:"),
      strong("Steffi Oesterreich (SO)
             Bahreini et al. Breast Cancer Research, 2017 #GSE89888"),
      helpText("Y537S/D538G mutations were knocked in MCF7 and T47D cells with rAAV and CRISPR/Cas9 respectively. 
               Genome-edited cells were hormone deprived and treated with or without 1 nM E2 for 24 hours. RNA were extraced and processed for RNA-seq."),
      strong("Simak Ali (SA)
             Harrod et al. Oncogene, 2017 #GSE78286"),
      helpText("Y537S mutation was knocked in MCF7 cells with CRISPR/Cas9. 
               Genome-edited cells were hormone deprived and treated with or without 1 nM E2 for 8 hours. RNA were extraced and processed for RNA-seq."),
      strong("Myles Brown (MB)
             Jeselsohn et al. Caner Cell, 2018 #GSE94493"),
      helpText("Dox-inducible Y537S or D538G plasmids were stably expressed in MCF7, Y537S mutant was expressed in T47D.
               Cells were hormone deprived and treated with or without 1 nM E2 for 24 hours. RNA were extraced and processed for RNA-seq."),         
               
    
    width=3),
    
    mainPanel(
      fluidRow(
        column(4,plotOutput("plotTPM_MCF7_SO")), 
        column(4,plotOutput("plotTPM_T47D_SO")), 
        column(4,plotOutput("plotTPM_MCF7_SA"))),
      fluidRow(
        column(4,plotOutput("plotTPM_T47D_MB")),
        column(4,plotOutput("plotTPM_MCF7_Y537S_MB")), 
        column(4,plotOutput("plotTPM_MCF7_D538G_MB")))
        
      ))
    
  
)
)