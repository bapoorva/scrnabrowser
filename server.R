library(shiny)
library(shinyBS)
library(RColorBrewer)
library(Biobase)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(shinyjs)
library(htmlwidgets)
library(DT)
library(shinyRGL)
library(rgl)
library(rglwidget)
library(Seurat)


server <- function(input, output,session) {
  
  ###################################################
  ###################################################
  ####### Display project list and load data  #######
  ###################################################
  ###################################################
  
  #Read the parameter file
  readexcel = reactive({
    file = read.csv("data/param.csv")
  })
  
  #Get Project list and populate drop-down
  output$projects = renderUI({
    excel=readexcel()
    prj=excel$projects
    selectInput("projects","Select a project",as.list(sort(as.character(prj))))
  })
  
  #display project list in Dashboard
  output$datasetTable<- renderTable({
    read.csv('data/param.csv',stringsAsFactors = F)
  }, digits = 1)
  
  #Load Rdata
  fileload <- reactive({
    inFile = paste('data/',as.character(input$projects),'.RData',sep = '')
    load(inFile)
    loaddata=scrna
    return(loaddata)
  })
  ###################################################
  ###################################################
  ####### Display Tsne plot with controls  ##########
  ###################################################
  ###################################################
  
  #get all categorical variables from metadata
  output$variables = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    metadata=metadata %>% select(starts_with("var"))
    var=colnames(metadata)
    selectInput("var","Select a Variable",var,"pick one")
  })
  
  tsneplot = reactive({
    scrna=fileload()
    variable=as.character(input$var)
    if(input$category=="clust"){
      TSNEPlot(object = scrna,group.by = "ident",do.hover = T,no.legend = FALSE)
    }else if(input$category=="var"){
      TSNEPlot(object = scrna,group.by = variable)
    }else if(input$category=="geneexp"){
      # validate(
      #   need(is.na(input$gene) == F, "Enter Gene name")
      # )
      genes=input$gene
      genes=unlist(strsplit(genes,","))
      FeaturePlot(object = scrna, features.plot = genes, cols.use = c("grey", "blue"),reduction.use = "tsne",
                  do.hover = T,no.legend = FALSE,data.hover = c("ident","nUMI", "nGene"))
    }
  })
  
  output$tsneplot = renderPlotly({
    input$category
    input$gene
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      tsneplot()
    })
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function(){
      paste0('tsneplot','.jpg',sep='')
    },
    content = function(file){
      png(file)
      jpeg(file, quality = 100, width = 800, height = 1300)
      tsneplot()
      dev.off()
    })
  
  ###################################################
  ###################################################
  ####### Display Biplot plot with controls #########
  ###################################################
  ###################################################

  output$bigeneplot <- renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    scrna=fileload()
    FeaturePlot(object = scrna, features.plot = c(input$bigene_genea,input$bigene_geneb), cols.use = c("grey","red","blue","green"),reduction.use = "tsne",
                no.legend = FALSE,overlay=TRUE,pt.size = input$bigene_pointsize)
    })
  })

  output$downloadbigene <- downloadHandler(
    filename = function(){
      paste0('biplot','.jpg',sep='')
    },
    content = function(file){
      png(file)
      jpeg(file, quality = 100, width = 800, height = 1300)
      bigeneplot()
      dev.off()
    })

  ###################################################
  ###################################################
  ####### Display Heatmap plot with controls  #######
  ###################################################
  ###################################################
  
  ###################################################
  ###################################################
  ####### Display DEG plot with controls  ###########
  ###################################################
  ###################################################
}#end of server