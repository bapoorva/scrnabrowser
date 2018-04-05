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
library(cowplot)


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
    metadata=as.data.frame(scrna@meta.data) %>% dplyr::select(-orig.ident,-old.ident,-res.0.6)
    #metadata=metadata %>% select(starts_with("var"))
    var=colnames(metadata)
    selectInput("var","Select a Variable",var,"pick one")
  })
  
  output$genes = renderUI({
    scrna=fileload()
    genelist=rownames(scrna@scale.data)
    selectizeInput("gene","Enter gene(s)",genelist, selected = "Hopx", multiple = FALSE,options = list(maxItems = 1))
  })
  
  tsneplot = reactive({
    scrna=fileload()
    pdf(NULL)
    metadata=as.data.frame(scrna@meta.data)
    metadata=metadata %>% select(starts_with("var"))
    variable=input$var
    feature=c("nGene","nUMI","percent.mito","S.Score","G2M.Score")
    tsne=c(colnames(metadata),"Phase")
    if(input$category =="clust"){
      TSNEPlot(object = scrna,group.by = "ident",do.hover = T,no.legend = FALSE,data.hover = c("ident","nUMI", "nGene"),do.label=TRUE)
    }else if(input$category=="geneexp"){
      genes=input$gene
      genes=unlist(strsplit(genes,","))
      FeaturePlot(object = scrna, features.plot = genes, cols.use = c("grey", "blue"),reduction.use = "tsne",
                  do.hover = T,no.legend = FALSE,data.hover = c("ident","nUMI", "nGene"))
    }else if(input$category =="var" & input$var %in% tsne){
      TSNEPlot(object = scrna,group.by = variable,do.hover = T,no.legend = FALSE,data.hover = variable)
    }else if(input$category =="var" & input$var %in% feature){
      FeaturePlot(object = scrna, features.plot = variable, cols.use = c("grey", "blue"),reduction.use = "tsne",do.hover = T,
                  no.legend = FALSE,data.hover = c("ident","nUMI", "nGene"))
    }
  })
  
  output$tsneplot = renderPlotly({
    input$category
    input$gene
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      tsneplot()
    })
  })
  
  tsneplot_dwld = reactive({
    scrna=fileload()
    pdf(NULL)
    metadata=as.data.frame(scrna@meta.data)
    metadata=metadata %>% select(starts_with("var"))
    variable=input$var
    feature=c("nGene","nUMI","percent.mito","S.Score","G2M.Score")
    tsne=c(colnames(metadata),"Phase")
    if(input$category =="clust"){
      TSNEPlot(object = scrna,group.by = "ident",no.legend = FALSE,do.label=TRUE,do.return=TRUE)
    }else if(input$category=="geneexp"){
      genes=input$gene
      genes=unlist(strsplit(genes,","))
      FeaturePlot(object = scrna, features.plot = genes, cols.use = c("grey", "blue"),reduction.use = "tsne",
                  no.legend = FALSE,do.return=TRUE)
    }else if(input$category =="var" & input$var %in% tsne){
      TSNEPlot(object = scrna,group.by = variable,no.legend = FALSE,do.return=TRUE)
    }else if(input$category =="var" & input$var %in% feature){
      FeaturePlot(object = scrna, features.plot = variable, cols.use = c("grey", "blue"),reduction.use = "tsne",
                  no.legend = FALSE,do.return=TRUE)
    }
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function(){
      paste0('tsneplot','.pdf',sep='')
    },
    content = function(file){
      pdf(file, width = 800, height = 1300)
      tsneplot_dwld()
      dev.off()
    })
  ###################################################
  ###################################################
  ####### Compare Tsne plots with controls  ##########
  ###################################################
  ###################################################
  output$tsnea = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data) %>% dplyr::select(-orig.ident,-old.ident,-res.0.6)
    #metadata=metadata %>% select(starts_with("var"))
    var=colnames(metadata)
    selectInput("tsnea","Select a Variable",var,"pick one")
  })
  
  output$tsneb = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data) %>% dplyr::select(-orig.ident,-old.ident,-res.0.6)
    #metadata=metadata %>% select(starts_with("var"))
    var=colnames(metadata)
    selectInput("tsneb","Select a Variable",var,"pick one")
  })
  
  comptsne = reactive({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    metadata=metadata %>% select(starts_with("var"))
    tsnea=input$tsnea
    tsneb=input$tsneb
    feature=c("nGene","nUMI","percent.mito","S.Score","G2M.Score")
    tsne=c(colnames(metadata),"Phase")
    if(input$categorya =="clust"){
      plot1=TSNEPlot(object = scrna,group.by = "ident",no.legend = FALSE,do.label = TRUE, do.return=T, pt.size = input$pointa)
    }else if(input$categorya=="geneexp"){
      # genes1=input$gene1
      # genes1=unlist(strsplit(genes1,","))
      plot1=FeaturePlot(object = scrna, features.plot = input$gene1, cols.use = c("grey", "blue"),reduction.use = "tsne",do.return=T,pt.size = input$pointa)
      plot1=eval(parse(text=paste("plot1$",input$gene1,sep="")))
    }else if(input$categorya =="var" & input$tsnea %in% tsne){
      plot1=TSNEPlot(object = scrna,group.by = tsnea,no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$pointa)
    }else if(input$categorya =="var" & input$tsnea %in% feature){
      plot1=FeaturePlot(object = scrna, features.plot = tsnea, cols.use = c("grey", "blue"),reduction.use = "tsne",do.return=T,pt.size = input$pointa)
      plot1=eval(parse(text=paste("plot1$",tsnea,sep="")))
    }

    if(input$categoryb =="clust"){
      plot2=TSNEPlot(object = scrna,group.by = "ident",no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$pointb)
    }else if(input$categoryb=="geneexp"){
      # genes2=input$gene2
      # genes2=unlist(strsplit(genes2,","))
      plot2=FeaturePlot(object = scrna, features.plot = input$gene2, cols.use = c("grey", "blue"),reduction.use = "tsne",do.return=T,pt.size = input$pointb)
      plot2=eval(parse(text=paste("plot2$",input$gene2,sep="")))
    }else if(input$categoryb =="var" & input$tsneb %in% tsne){
      plot2=TSNEPlot(object = scrna,group.by = tsneb,no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$pointb)
    }else if(input$categoryb =="var" & input$tsneb %in% feature){
      plot2=FeaturePlot(object = scrna, features.plot = tsneb, cols.use = c("grey", "blue"),reduction.use = "tsne",do.return=T,pt.size = input$pointb)
      plot2=eval(parse(text=paste("plot2$",tsneb,sep="")))
    }

    plot_grid(plot1,plot2)
    
    # tsnea<- TSNEPlot(object = scrna,do.label = TRUE, do.return=T,label.size = 3)
    # tsneb<- FeaturePlot(object = scrna, features.plot ="nUMI", cols.use = c("grey", "blue"),reduction.use = "tsne",do.return=T)
    # plot_grid(tsnea,tsneb$nUMI)

  })
  
  output$comptsne = renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      comptsne()
    })
  })
  ###################################################
  ###################################################
  ####### Display Biplot plot with controls #########
  ###################################################
  ###################################################
  
  bigeneplot <- reactive({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      scrna=fileload()
      FeaturePlot(object = scrna, features.plot = c(input$bigene_genea,input$bigene_geneb), cols.use = c("grey","red","blue","green"),reduction.use = "tsne",
                  no.legend = FALSE,overlay=TRUE,pt.size = input$bigene_pointsize,do.return = T)
    })
  })
  
  output$bigeneplot <- renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    bigeneplot()
    })
  })

  output$downloadbigene <- downloadHandler(
    filename = function(){
      paste0('biplot','.jpg',sep='')
    },
    content = function(file){
      jpeg(file, quality = 100, width = 800, height = 1300)
      plot(bigeneplot())
      dev.off()
    })
  
  ###################################################
  ###################################################
  ####### Display DEG plot with controls  ###########
  ###################################################
  ###################################################
  markergenes = reactive({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    scrna=fileload()
    if(input$identb==""){
      markers=FindMarkers(object = scrna, ident.1 = input$identa, min.pct = input$minpct,logfc.threshold=input$lfc,test.use=input$test)
    }else{
      identb=input$identb
      p=unlist(strsplit(identb,","))
    markers=FindMarkers(object = scrna, ident.1 = input$identa, ident.2 = p, min.pct = input$minpct,logfc.threshold=input$lfc,test.use=input$test)
    }
    })
  })
  
  output$markergenes = DT::renderDataTable({
    input$identa
    input$identb
    DT::datatable(markergenes(),
                  extensions = 'Buttons', options = list(
                    dom = 'Bfrtip',
                    buttons = list()),
                  rownames=TRUE,selection = list(mode = 'single', selected =1),escape=FALSE)
  })
  
  ###################################################
  ###################################################
  ####### Display feature plot with controls  #######
  ###################################################
  ###################################################
  
  output$markergenes_out = DT::renderDataTable({
    input$identa
    input$identb
    DT::datatable(markergenes(),
                  extensions = 'Buttons', options = list(
                    dom = 'Bfrtip',
                    buttons = list()),
                  rownames=TRUE,selection = list(mode = 'single', selected =1),escape=FALSE)
  })
  
  markgeneplot <- reactive({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      scrna=fileload()
      markers=markergenes()
      s=input$markergenes_out_rows_selected # get  index of selected row from table
      markers=markers[s, ,drop=FALSE]
      FeaturePlot(object = scrna, features.plot = rownames(markers), cols.use = c("grey","blue"),reduction.use = "tsne",
                  no.legend = FALSE,pt.size = input$marker_pointsize,do.return = T)
    })
  })
  
  # markgeneplot <- reactive({
  #   withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
  #     scrna=fileload()
  #     markers=markergenes()
  #     FeaturePlot(object = scrna, features.plot = rownames(markers)[1:input$cowplot], cols.use = c("grey","blue"),reduction.use = "tsne",
  #                 no.legend = FALSE,pt.size = input$marker_pointsize)
  #   })
  # })
  
  output$markgeneplot <- renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      markgeneplot()
    })
  })
  
  output$downloaddeg <- downloadHandler(
    filename = function() { paste(input$projects, '.csv', sep='') },
    content = function(file) {
      write.csv(markergenes(), file)
    })
  
  output$dwldfeature <- downloadHandler(
    filename = function(){
      paste0('featureplot','.jpg',sep='')
    },
    content = function(file){
      jpeg(file, quality = 100, width = 800, height = 1300)
      markgeneplot()
      dev.off()
    })
  
  ###################################################
  ###################################################
  ####### Display Violin  plot with controls  #######
  ###################################################
  ###################################################
  output$grptype = renderUI({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    metadata=metadata %>% select(starts_with("var"))
    var=c("ident",colnames(metadata))
    selectInput("grptype","Select a Variable",var,"pick one")
    })
  })
  
   violinplot <- reactive({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      scrna=fileload()
      markers=markergenes()
      VlnPlot(object = scrna, features.plot = rownames(markers)[1:input$vplot],group.by = input$grptype,do.return = T)
    })
  })
  
   output$violinplot <- renderPlot({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       violinplot()
     })
   })
  ###################################################
  ###################################################
  ####### Display Heatmap plot with controls  #######
  ###################################################
  ###################################################
   output$heatmapgenes = renderUI({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       markers=markergenes()
       validate(
         need(nrow(markers)>0, "No Marker genes found")
       )
       if(nrow(markers)<10){
         min=1
         max=nrow(markers)
       }else{
         min=10
         max=nrow(markers)
       }
       sliderInput("heatmapgenes", "Number of top genes to plot:",min = min, max = max,value = min)
     })
   })
   
   output$hmpgrp = renderUI({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       scrna=fileload()
       metadata=as.data.frame(scrna@meta.data)
       metadata=metadata %>% select(starts_with("var"))
       var=c("ident",colnames(metadata))
       selectInput("hmpgrp","Select a Variable",var,"pick one")
     })
   })
   
   output$heatmap <- renderPlot({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       scrna=fileload()
       markers=markergenes()
       if(input$hmpcol=="PuYl"){
         lowcol="darkmagenta"
         midcol="black"
         highcol="yellow"
       }else if(input$hmpcol=="BuGn"){
         lowcol="yellow"
           midcol="green"
           highcol="blue"
       }else if(input$hmpcol=="RdYl"){
         lowcol="yellow"
           midcol="red"
           highcol="black"
       }else if(input$hmpcol=="RdBu"){
         lowcol="red"
         midcol="white"
         highcol="blue"}
       DoHeatmap(object = scrna, genes.use = rownames(markers)[1:input$heatmapgenes],group.by = input$hmpgrp, draw.line= T,
                 group.label.rot= T, col.low=lowcol, col.mid =midcol ,col.high = highcol,slim.col.label=TRUE)
     })
   })
  
}#end of server