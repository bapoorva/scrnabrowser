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
    prj=as.character(excel$projects)
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
  output$tsnea2 = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data) 
    #metadata=metadata %>% select(starts_with("var"))
    var=colnames(metadata)
    selectInput("tsnea2","Select a Variable",var,"pick one")
  })
  
  output$tsneb2 = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data) 
    #metadata=metadata %>% select(starts_with("var"))
    var=colnames(metadata)
    selectInput("tsneb2","Select a Variable",var,"pick one")
  })
  
  output$subsaui = renderUI({
    scrna=fileload()
    clusts=levels(scrna@ident)
    if(input$categorya2=="clust"){
      selectInput("selclust","Select a Cluster",clusts)
    }else if(input$categorya2=="var"){
      metadata=as.data.frame(scrna@meta.data)
      met= sapply(metadata,is.numeric)
      feature=names(met[met==TRUE])
      tsne=names(met[met==FALSE])
      t=paste("scrna@meta.data$",input$tsnea2,sep="")
      if(input$tsnea2 %in% tsne){
        opt1=unique(eval(parse(text=t)))
        selectInput("selclust2","Select one of the options",opt1)
      }else if(input$tsnea2 %in% feature){
        min=min(eval(parse(text=t)))
        max=max(eval(parse(text=t)))
        sliderInput("tsnea2lim", label = h5("Select Range"), min = min,max =max, value =c(min,max)) 
      }
    }else if(input$categorya2=="geneexp"){
      validate(need(input$categorya2!="geneexp","Cannot subselect gene expression values"))
    }
  })
  
  output$subsbui = renderUI({
    scrna=fileload()
    clusts=levels(scrna@ident)
    if(input$categoryb2=="clust"){
      selectInput("selclustb","Select a Cluster",clusts)
    }else if(input$categoryb2=="var"){
      metadata=as.data.frame(scrna@meta.data)
      met= sapply(metadata,is.numeric)
      feature=names(met[met==TRUE])
      tsne=names(met[met==FALSE])
      t=paste("scrna@meta.data$",input$tsneb2,sep="")
      if(input$tsneb2 %in% tsne){
        opt2=unique(eval(parse(text=t)))
        selectInput("selclustb2","Select one of the options",opt2)
      }else if(input$tsneb2 %in% feature){
        min=min(eval(parse(text=t)))
        max=max(eval(parse(text=t)))
        sliderInput("tsneb2lim", label = h5("Select Range"), min = min,max =max, value =c(min,max)) 
      }
    }else if(input$categoryb2=="geneexp"){
      validate(need(input$categoryb2!="geneexp","Cannot subselect gene expression values"))
    }
  })

  
  comptsne2 = reactive({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    met= sapply(metadata,is.numeric)
    #metadata=metadata %>% select(starts_with("var"))
    tsnea=input$tsnea2
    tsneb=input$tsneb2
    feature=names(met[met==TRUE])
    #feature=c("nGene","nUMI","percent.mito","S.Score","G2M.Score","var.ratio.pca")
    tsne=names(met[met==FALSE])
    #tsne=c(colnames(metadata),"Phase","sample")
    if(input$categorya2 =="clust" & input$subsa==F){
      plot1=TSNEPlot(object = scrna,group.by = "ident",no.legend = FALSE,do.label = TRUE, do.return=T, pt.size = input$pointa2)
    }else if(input$categorya2 =="clust" & input$subsa==TRUE){
      cells=names(scrna@ident[scrna@ident==input$selclust])
      plot1=TSNEPlot(object = scrna,cells.highlight=cells,group.by = "ident",no.legend = FALSE,do.label = TRUE, do.return=T, pt.size = input$pointa2)
    }else if(input$categorya2=="geneexp"){
      plot1=FeaturePlot(object = scrna, features.plot = input$gene1a, cols.use = c("grey", "blue"),reduction.use = "tsne",do.return=T,pt.size = input$pointa2,no.legend = FALSE)
      plot1=eval(parse(text=paste("plot1$",input$gene1a,sep="")))
    }else if(input$categorya2 =="var" & input$tsnea2 %in% tsne & input$subsa==FALSE){
      plot1=TSNEPlot(object = scrna,group.by = tsnea,no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$pointa2)
    }else if(input$categorya2 =="var" & input$tsnea2 %in% tsne & input$subsa==TRUE){
      t=paste("rownames(scrna@meta.data[scrna@meta.data$",input$tsnea2,"==\"",input$selclust2,"\",])",sep="")
      cells=eval(parse(text=t))
      plot1=TSNEPlot(object = scrna,group.by = tsnea,cells.highlight=cells,no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$pointa2)
    }else if(input$categorya2 =="var" & input$tsnea2 %in% feature & input$subsa==FALSE){
      plot1=FeaturePlot(object = scrna, features.plot = tsnea, cols.use = c("grey", "blue"),reduction.use = "tsne",do.return=T,pt.size = input$pointa2,no.legend = FALSE)
      plot1=eval(parse(text=paste("plot1$",tsnea,sep="")))
    }else if(input$categorya2 =="var" & input$tsnea2 %in% feature & input$subsa==TRUE){
      t=paste('rownames(scrna@meta.data[scrna@meta.data$',input$tsnea2, '>',input$tsnea2lim[1], ' & metadata$',input$tsnea2, '<', input$tsnea2lim[2],',])',sep="")
      cells=eval(parse(text=t))
      plot1=FeaturePlot(object = scrna, features.plot = tsnea,cells.use = cells, cols.use = c("grey", "blue"),reduction.use = "tsne",do.return=T,pt.size = input$pointa2,no.legend = FALSE)
      plot1=eval(parse(text=paste("plot1$",tsnea,sep="")))
    }
    
    if(input$categoryb2 =="clust" & input$subsb==F){
      plot2=TSNEPlot(object = scrna,group.by = "ident",no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$pointa2)
    }else if(input$categoryb2 =="clust" & input$subsb==TRUE){
      cells=names(scrna@ident[scrna@ident==input$selclustb])
      plot2=TSNEPlot(object = scrna,cells.highlight=cells,group.by = "ident",no.legend = FALSE,do.label = TRUE, do.return=T, pt.size = input$pointa2)
    }else if(input$categoryb2=="geneexp"){
      plot2=FeaturePlot(object = scrna, features.plot = input$gene2a, cols.use = c("grey", "blue"),reduction.use = "tsne",do.return=T,pt.size = input$pointa2,no.legend = FALSE)
      plot2=eval(parse(text=paste("plot2$",input$gene2a,sep="")))
    }else if(input$categoryb2 =="var" & input$tsneb2 %in% tsne & input$subsb==F){
      plot2=TSNEPlot(object = scrna,group.by = tsneb,no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$pointa2)
    }else if(input$categoryb2 =="var" & input$tsneb2 %in% tsne & input$subsb==TRUE){
      t=paste("rownames(scrna@meta.data[scrna@meta.data$",input$tsneb2,"==\"",input$selclustb2,"\",])",sep="")
      cells=eval(parse(text=t))
      plot2=TSNEPlot(object = scrna,group.by = tsneb,cells.highlight=cells,no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$pointa2)
    }else if(input$categoryb2 =="var" & input$tsneb2 %in% feature & input$subsb==F){
      plot2=FeaturePlot(object = scrna, features.plot = tsneb, cols.use = c("grey", "blue"),reduction.use = "tsne",do.return=T,pt.size = input$pointa2,no.legend = FALSE)
      plot2=eval(parse(text=paste("plot2$",tsneb,sep="")))
    }else if(input$categoryb2 =="var" & input$tsneb2 %in% feature & input$subsb==TRUE){
      t=paste('rownames(scrna@meta.data[scrna@meta.data$',input$tsneb2, '>',input$tsneb2lim[1], ' & metadata$',input$tsneb2, '<', input$tsneb2lim[2],',])',sep="")
      cells=eval(parse(text=t))
      plot2=FeaturePlot(object = scrna, features.plot = tsneb,cells.use = cells, cols.use = c("grey", "blue"),reduction.use = "tsne",do.return=T,pt.size = input$pointa2,no.legend = FALSE)
      plot2=eval(parse(text=paste("plot2$",tsneb,sep="")))
    }
    
    plot_grid(plot1,plot2)

  })
  
  output$comptsne2 = renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      comptsne2()
    })
  })
  
  output$downloadtsneplot <- downloadHandler(
    filename = function() {
      paste0("Compare_tsne.jpg")
    },
    content = function(file){
      ggsave(file, plot = comptsne2(), device = "jpg")
    })
  ###################################################
  ###################################################
  ####### Compare Tsne plots with controls  ##########
  ###################################################
  ###################################################
  output$tsnea = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data) 
    #metadata=metadata %>% select(starts_with("var"))
    var=colnames(metadata)
    selectInput("tsnea","Select a Variable",var,"pick one")
  })
  
  output$tsneb = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data) 
    #metadata=metadata %>% select(starts_with("var"))
    var=colnames(metadata)
    selectInput("tsneb","Select a Variable",var,"pick one")
  })
  
  comptsne = reactive({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    met= sapply(metadata,is.numeric)
    #metadata=metadata %>% select(starts_with("var"))
    tsnea=input$tsnea
    tsneb=input$tsneb
    feature=names(met[met==TRUE])
    #feature=c("nGene","nUMI","percent.mito","S.Score","G2M.Score","var.ratio.pca")
    tsne=names(met[met==FALSE])
    #tsne=c(colnames(metadata),"Phase","sample")
    if(input$categorya =="clust"){
      plot1=TSNEPlot(object = scrna,group.by = "ident",no.legend = FALSE,do.label = TRUE, do.return=T, pt.size = input$pointa)
    }else if(input$categorya =="var" & input$tsnea %in% tsne){
      plot1=TSNEPlot(object = scrna,group.by = tsnea,no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$pointa)
    }else if(input$categorya =="var" & input$tsnea %in% feature){
      plot1=FeaturePlot(object = scrna, features.plot = tsnea, cols.use = c("grey", "blue"),reduction.use = "tsne",do.return=T,pt.size = input$pointa)
      plot1=eval(parse(text=paste("plot1$",tsnea,sep="")))
    }

    markers=markergenes()
      s=input$markergenes_rows_selected # get  index of selected row from table
      markers=markers[s, ,drop=FALSE]
      plot2=FeaturePlot(object = scrna, features.plot = rownames(markers), cols.use = c("grey","blue"),reduction.use = "tsne",
                        no.legend = FALSE,pt.size = input$pointa,do.return = T)
      plot2=eval(parse(text=paste("plot2$",rownames(markers),sep="")))
    
      plot3=VlnPlot(object = scrna, features.plot = rownames(markers),group.by = input$grptype,do.return = T)
      plot4=RidgePlot(object = scrna, features.plot = rownames(markers),group.by = input$grptype,do.return = T)
      
    
      row1=plot_grid(plot1,plot2,align = 'h', rel_heights = c(1, 1),axis="lr", nrow=1)
      row2=plot_grid(plot3,plot4,align = 'h', rel_heights = c(1, 1),axis="lr", nrow=1)
    plot_grid(row1,row2,align = 'v', rel_heights = c(1.7, 1),axis="tb",ncol=1)

  })
  
  output$comptsne = renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      comptsne()
    })
  })
  
  output$downloadplot <- downloadHandler(
    filename = function() {
      paste0("Plot.jpg")
    },
    content = function(file){
      ggsave(file, plot = comptsne(), device = "jpg", height = 1000 ,width =700 )
    })
  ###################################################
  ###################################################
  ####### Display Biplot plot with controls #########
  ###################################################
  ###################################################
  
  # bigeneplot <- reactive({
  #   withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
  #     scrna=fileload()
  #     FeaturePlot(object = scrna, features.plot = c(input$bigene_genea,input$bigene_geneb), cols.use = c("grey","red","blue","green"),reduction.use = "tsne",
  #                 no.legend = FALSE,overlay=TRUE,pt.size = input$bigene_pointsize,do.return = T)
  #   })
  # })
  # 
  # output$bigeneplot <- renderPlot({
  #   withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
  #   bigeneplot()
  #   })
  # })
  # 
  # output$downloadbigene <- downloadHandler(
  #   filename = function(){
  #     paste0('biplot','.jpg',sep='')
  #   },
  #   content = function(file){
  #     jpeg(file, quality = 100, width = 800, height = 1300)
  #     plot(bigeneplot())
  #     dev.off()
  #   })
  ######################################################################################################
  ######################################################################################################
  ####### Display Biplot plot with controls ############################################################
  ######################################################################################################
  ######################################################################################################
  getGeneRange <- function(scrna,gene_probes){
    gene_values=as.data.frame(FetchData(scrna,gene_probes[1]))
    minr<- round(min(gene_values),2) 
    maxr<- round(max(gene_values),2)
    
    return(c(ifelse(minr==0,.1,minr-.1),maxr))
    #return(c(1,12))
    
  }
  
  output$bigene_rangea <- renderUI({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    r<-getGeneRange(fileload(),input$bigene_genea)
    sliderInput("bigene_rangea", "Expression Limit Gene A (log2 UMI)",
                min = 0, max = r[2], value = c(r[1],r[2]),step=.25)
  })
  })
  
  output$bigene_rangeb <- renderUI({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    r<-getGeneRange(fileload(),input$bigene_geneb)
    sliderInput("bigene_rangeb", "Expression Limit Gene B (log2 UMI)",
                min = 0, max = r[2], value = c(r[1],r[2]),step=.25)
  })
  })
  
  output$bigeneplot <- renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    bigene_plot(fileload(),
                c(input$bigene_genea,input$bigene_geneb),
                limita=input$bigene_rangea,
                limitb=input$bigene_rangeb,
                marker_size = input$bigene_pointsize)
  })
  })

  bigene_getValues <- function(scrna,gene_probes,limita,limitb){
    gene_values=FetchData(scrna,c(gene_probes[1],gene_probes[2]))
    colnames(gene_values) <- c('genea','geneb')
    # gene_values=as.data.frame(gene_values) %>% 
    #   mutate(value =ifelse(genea>=limita[1] & geneb <limitb[1], gene_probes[1], ifelse(genea<limita[1] & geneb >=limitb[1],
    #                   gene_probes[2],ifelse(genea>=limita[1] & geneb >=limitb[1],"DoublePos","NULL"))))
    as.data.frame(gene_values) %>% 
      mutate(value = ifelse(genea>=limita[1] & geneb>=limitb[1],
                            'both',
                            ifelse(genea>=limita[1] & geneb<limitb[1],
                                   gene_probes[1],
                                   ifelse(genea<=limita[1] & geneb>=limitb[1],
                                          gene_probes[2],
                                          'none')))
      ) #%>% select(value)
  }
  
  monocle_theme_opts <- function()
  {
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
      theme(panel.border = element_blank()) +
      theme(axis.line.x = element_line(size=0.25, color="black")) +
      theme(axis.line.y = element_line(size=0.25, color="black")) +
      theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
      theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
      theme(panel.background = element_rect(fill='white')) +
      theme(legend.key=element_blank())
  }
  
  bigene_plot <- function (scrna, gene_probes, x=1,y=2, limita=c(1,100), limitb=c(1,100), marker_size = 0.1,
                           title = NULL)
  {

    gene_values <- bigene_getValues(scrna,gene_probes,limita,limitb)
    projection=as.data.frame(scrna@dr$tsne@cell.embeddings)
    colnames(projection) <- c("Component.1", "Component.2")
    proj_gene <- data.frame(cbind(projection, gene_values))
    #proj_gene$value = factor(proj_gene$value,levels=unique(proj_gene$value))
    proj_gene$value = factor(proj_gene$value,levels=c('both',gene_probes[1],gene_probes[2],'none'))
    proj_gene <- arrange(proj_gene, desc(value))

    p <- ggplot(proj_gene, aes(Component.1, Component.2)) +
      geom_point(aes(colour = value), size = marker_size) +
      scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A", 'grey90'),drop=F) +
      theme(legend.key.size = unit(10,"point")) + xlab(paste("Component", x)) +
      ylab(paste("Component", y))


    if (!is.null(title)) {
      p <- p + ggtitle(title)
    }
    p <- p + monocle_theme_opts() + theme(plot.title = element_text(hjust = 0.5),
                                          legend.position="bottom",
                                          legend.title=element_blank(),
                                          legend.text=element_text(size=14),
                                          panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank())
    return(p)
  }

  ####################################################
  ###################################################
  ########## Set ident to choose markers  ###########
  ###################################################
  ###################################################
  output$setidentlist = renderUI({
      scrna=fileload()
      metadata=as.data.frame(scrna@meta.data)
      metadata=metadata %>% select(starts_with("var_"))
      var=c(colnames(metadata))
      selectInput("setidentlist","Choose category to compare",var,"pick one")
    
  })
  
  output$identa = renderUI({
    scrna=fileload()
    if(input$setident==T){
      scrna <- SetAllIdent(object = scrna, id = input$setidentlist)
      options=unique(scrna@ident)
    }else{
      options=as.numeric(levels(scrna@ident))
    }
    selectInput("identa", "First cluster to compare",options)
    })
  
    output$identb = renderUI({
      scrna=fileload()
      if(input$setident==T){
        scrna <- SetAllIdent(object = scrna, id = input$setidentlist)
        options=unique(scrna@ident)
      }else{
        options=as.numeric(levels(scrna@ident))
      }
        selectInput("identb", "Second cluster to compare",options,selected=options[2])
    })
  
  
  ###################################################
  ###################################################
  ####### Display DEG plot with controls  ###########
  ###################################################
  ###################################################
  markergenes = reactive({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    scrna=fileload()
    if(input$setident==T){
      scrna <- SetAllIdent(object = scrna, id = input$setidentlist)
    }
    
    if(input$identb==""){
      markers=FindMarkers(object = scrna, ident.1 = input$identa, min.pct = input$minpct,logfc.threshold=input$lfc,test.use=input$test)
      geneid=rownames(markers)
      url= paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",geneid,sep = "")
      markers$Link=paste0("<a href='",url,"'target='_blank'>",rownames(markers),"</a>")
      
    }else{
      identb=input$identb
      p=unlist(strsplit(identb,","))
    markers=FindMarkers(object = scrna, ident.1 = input$identa, ident.2 = p, min.pct = input$minpct,logfc.threshold=input$lfc,test.use=input$test)
    geneid=rownames(markers)
    url= paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",geneid,sep = "")
    markers$Link=paste0("<a href='",url,"'target='_blank'>",rownames(markers),"</a>")
    
    }
    
    })
    return(markers)
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
  
  output$downloaddeg <- downloadHandler(
    filename = function() { paste(input$projects, '.csv', sep='') },
    content = function(file) {
      write.csv(markergenes(), file)
    })
  

  
  ####################################################
  ###################################################
  ####### Display Violin  plot with controls  #######
  ###################################################
  ###################################################
  output$grptype = renderUI({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    metadata=metadata %>% select(starts_with("var_"))
    var=c("ident",colnames(metadata))
    selectInput("grptype","Select a Variable to group the Violin Plot",var,"pick one")
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
   
   output$downloadheatmap <- downloadHandler(
     filename = function() {
       paste0("heatmap.jpg")
     },
     content = function(file){
       jpeg(file, quality = 100, width = 800, height = 400)
       plot(heatmap())
       dev.off()
     })
  
   ###################################################
   ###################################################
   ####### UPLOAD GENELISTS AND SELECT CLUSTER #######
   ###################################################
   ###################################################
   output$list1 <- renderUI({
     fileInput('genelist1', 'Upload Receptor Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
   })
   
   output$list2 <- renderUI({
     fileInput('genelist2', 'Upload Ligand Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
   })
   
   output$clust1 <- renderUI({
     selectInput("clust1","Pick cluster1",c(0:11),selected=1)
   })
   
   output$clust2 <- renderUI({
     selectInput("clust2","Pick cluster2",c(0:11),selected=2)
   })
   
   output$list1.1 <- renderUI({
     fileInput('genelist1.1', 'Upload Gene List1',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
   })
   
   output$list2.1 <- renderUI({
     fileInput('genelist2.1', 'Upload Gene List2',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
   })
   
   output$clust1.1 <- renderUI({
     selectInput("clust1.1","Pick cluster1",c(0:11),selected=1)
   })
   
   output$clust2.1 <- renderUI({
     selectInput("clust2.1","Pick cluster2",c(0:11),selected=2)
   })
   ###################################################
   ###################################################
   #### Load lig-receptor list and display results ###
   ###################################################
   ###################################################
   
   datasetInput = reactive({
     scrna=fileload()
     tt=rownames(scrna@raw.data)
     genes=fread("data/ligrecgenes.txt",header = TRUE)
     genes2=tt[tt %in% genes$genes]
     my.data=FetchData(scrna,c("ident","nGene",genes2))
     my.data= my.data %>% rename('clust'='ident')
     result=data.frame()
     res=data.frame()
     for(i in 0:(length(unique(my.data$clust))-1)){
       for(j in 0:(length(unique(my.data$clust))-1)){
         if(i!=j){
           test=my.data[my.data$clust==i | my.data$clust==j,]
           R_c1=test[test$clust==i ,(colnames(test) %in% rl$Receptor.ApprovedSymbol)]
           L_c2=test[test$clust==j , (colnames(test) %in% rl$Ligand.ApprovedSymbol)]
           keep1 = colSums(R_c1>1)>=.5*dim(R_c1)[1]
           keep2 = colSums(L_c2>1)>=.5*dim(L_c2)[1]
           
           R_c1=R_c1[,keep1]
           L_c2=L_c2[,keep2]
           res=rl[(rl$Ligand.ApprovedSymbol %in% colnames(L_c2)) & (rl$Receptor.ApprovedSymbol %in% colnames(R_c1)),]
           
         }
         else{}
         if(nrow(res)!=0){
           res$Receptor_cluster=i
           res$Lig_cluster=j
           result=rbind(result,res)
         }else{result=result}
       }
     }
     result=result[result$Receptor_cluster!=result$Lig_cluster,]
     
     if(input$clust=="clust" & input$gene=="allgene"){
       clusters=c(input$clust1,input$clust2)
       result=result[(result$Receptor_cluster %in% clusters) & (result$Lig_cluster%in% clusters),]
     }else if(input$clust=="clust" & input$gene=="genelist"){
       clusters.1=c(input$clust1.1,input$clust2.1)
       result=result[(result$Receptor_cluster %in% clusters.1) & (result$Lig_cluster%in% clusters.1),]
     }else{result=result
     }
     
     if(input$gene=="genelist"){
       if(input$clust=="all"){
         g1=input$genelist1
         g2=input$genelist2
       }else if(input$clust=="clust"){
         g1=input$genelist1.1
         g2=input$genelist2.1
       }
       genes=read.table(g1$datapath,stringsAsFactors = F)#get complete gene list as string
       g1=as.vector(genes$V1)
       g1=tolower(g1)
       firstup <- function(x) {
         substr(x, 1, 1) <- toupper(substr(x, 1, 1))
         x
       }
       g1=firstup(g1)
       genes2=read.table(g2$datapath,stringsAsFactors = F)
       g2=as.vector(genes2$V1)
       g2=tolower(g2)
       g2=firstup(g2)
       result=result[(result$Receptor.ApprovedSymbol %in% g1) & (result$Ligand.ApprovedSymbol %in% g2),]
     }else{
       result=result
     }
     return(result)
   })
   
   #print data TABLE
   output$pairs_res = DT::renderDataTable({
     input$clust1
     input$clust2
     input$genelist1
     input$genelist2
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
       DT::datatable(datasetInput(),
                     extensions = c('Buttons','Scroller'),
                     options = list(dom = 'Bfrtip',
                                    searchHighlight = TRUE,
                                    pageLength = 10,
                                    lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                    scrollX = TRUE,
                                    buttons = c('copy', 'print')
                     ),rownames=FALSE,caption= "Result")
     })
   })
   
   
   ###################################################
   ###################################################
   ####### LOAD EXCEL AND POPULATE DROP DOWN #########
   ###################################################
   ###################################################
   
   #Read the parameter file
   readexcel = reactive({
     file = read.csv("data/param.csv")
   })
   
   #Get Project list and populate drop-down
   output$projects = renderUI({
     excel=readexcel()
     excel=excel[excel$type=="scrna",]
     prj=excel$projects
     selectInput("projects","Select a project",as.list(sort(as.character(prj))))
   })
   
   output$recprj = renderUI({
     excel=readexcel()
     excel=excel[excel$type==input$receptor,]
     prj=excel$projects
     selectInput("recprj","Select a project",as.list(sort(as.character(prj))))
   })
   
   output$ligprj = renderUI({
     excel=readexcel()
     excel=excel[excel$type==input$ligand,]
     prj=excel$projects
     selectInput("ligprj","Select a project",as.list(sort(as.character(prj))))
   })
   ###################################################
   ###################################################
   ############# READ OR LOAD DATA FILES #############
   ###################################################
   ###################################################
   recfileload <- reactive({
     inFile = paste('data/',as.character(input$recprj),'.RData',sep = '')
     load(inFile)
     loaddata=results
     return(loaddata)
   })
   ligfileload <- reactive({
     inFile = paste('data/',as.character(input$ligprj),'.RData',sep = '')
     load(inFile)
     loaddata=results
     return(loaddata)
   })
   
   recfileread <- reactive({
     inFile = paste('data/',as.character(input$recprj),'.csv',sep = '')
     loaddata=read.csv(inFile)
     return(loaddata)
   })
   
   ligfileread <- reactive({
     inFile = paste('data/',as.character(input$ligprj),'.csv',sep = '')
     loaddata=read.csv(inFile)
     return(loaddata)
   })
   
   ###################################################
   ###################################################
   ### POPULATE DROPDOWN FOR MAINEFFECT AND CLUSTER ##
   ###################################################
   ###################################################
   output$rectype = renderUI({
     if(input$receptor=="scrna"){
       loaddata=recfileread()
       maineffect=unique(loaddata$ident)
     }else if(input$receptor=="rna" | input$receptor=="microarray"){
       results=recfileload()
       pd=pData(results$eset)
       maineffect=pd$maineffect
     }
     selectInput("rectype","Select subtype (maineffect for RNAseq and cluster for scRNA)",as.list(sort(as.character(maineffect))),selected = 1)
   })
   
   output$ligtype = renderUI({
     if(input$ligand=="scrna"){
       loaddata=ligfileread()
       maineffect=unique(loaddata$ident)
     }else if(input$ligand=="rna" | input$ligand=="microarray"){
       results=ligfileload()
       pd=pData(results$eset)
       maineffect=pd$maineffect
     }
     selectInput("ligtype","Select subtype (maineffect for RNAseq and cluster for scRNA)",as.list(sort(as.character(maineffect))),selected = 1)
   })
   
   ###################################################
   ###################################################
   ####### COMPARE EFFECTS AND DISPLAY RESULTS #######
   ###################################################
   ###################################################
   firstup <- function(x) {
     substr(x, 1, 1) <- toupper(substr(x, 1, 1))
     x
   }
   
   ligrecpairs = reactive({
     if(input$receptor=="scrna"){
       scrna=recfileread()
       scrna=scrna[scrna$ident==input$rectype,]
       rownames(scrna)=scrna$X
       scrna=scrna %>% dplyr::select(-X:-nGene)
       recumi=input$recumi
       recsamp=(input$recsamp)/100
       keep= colSums(scrna>recumi)>=recsamp*dim(scrna)[1]
       scrna2=scrna[,keep]
       rec_avg=NULL
       rec_genes=colnames(scrna2)
     }else if(input$receptor=="rna" | input$receptor=="microarray"){
       results=recfileload()
       loaddata=as.data.frame(exprs(results$eset))
       pd=pData(results$eset)
       #minexpr=unique(pd$minexpr)
       minexpr=input$exprec
       sel=input$rectype
       samp=as.character(pd$sample_name[pd$maineffect==sel])
       loaddata2=loaddata %>% dplyr::select(samp)
       loaddata2$avg=rowMeans(loaddata2)
       fd=fData(results$eset)
       fd$id=rownames(fd)
       loaddata2$id=rownames(loaddata2)
       loaddata2=inner_join(loaddata2,fd,by="id")
       loaddata2=loaddata2[loaddata2$avg > minexpr,]
       rec_avg=loaddata2 %>% dplyr::select(avg,SYMBOL)
       rec_genes=loaddata2$SYMBOL
     }
     
     if(input$ligand=="scrna"){
       scrna=ligfileread()
       scrna=scrna[scrna$ident==input$ligtype,]
       rownames(scrna)=scrna$X
       scrna=scrna %>% dplyr::select(-X:-nGene)
       ligumi=input$ligumi
       ligsamp=(input$ligsamp)/100
       keep= colSums(scrna>ligumi)>=ligsamp*dim(scrna)[1]
       scrna2=scrna[,keep]
       lig_avg=NULL
       lig_genes=colnames(scrna2)
     }else if(input$ligand=="rna" | input$ligand=="microarray"){
       results=ligfileload()
       loaddata=as.data.frame(exprs(results$eset))
       pd=pData(results$eset)
       #minexpr=unique(pd$minexpr)
       minexpr=input$explig
       sel=input$ligtype
       samp=as.character(pd$sample_name[pd$maineffect==sel])
       loaddata2=loaddata %>% dplyr::select(samp)
       loaddata2$avg=rowMeans(loaddata2)
       fd=fData(results$eset)
       fd$id=rownames(fd)
       loaddata2$id=rownames(loaddata2)
       loaddata2=inner_join(loaddata2,fd,by="id")
       loaddata2=loaddata2[loaddata2$avg > minexpr,]
       lig_avg=loaddata2 %>% dplyr::select(avg,SYMBOL)
       lig_genes=loaddata2$SYMBOL
     }
     rl=read.csv("data/lig-rec.csv")
     if(is.null(lig_avg)==F){rl=left_join(rl,lig_avg,by=c("Ligand.ApprovedSymbol"="SYMBOL")) %>% dplyr::rename(Ligand_AvgExpr=avg)}
     if(is.null(rec_avg)==F){rl=left_join(rl,rec_avg,by=c("Receptor.ApprovedSymbol"="SYMBOL")) %>% dplyr::rename(Receptor_AvgExpr=avg)}
     rl=rl[rl$Ligand.ApprovedSymbol %in% lig_genes & rl$Receptor.ApprovedSymbol %in% rec_genes,]
     
     
     if(input$liggene ==T){
       validate(
         need(input$liggeneli, "Please Upload genelist")
       )
       lgene=input$liggeneli
       lgenes=read.table(lgene$datapath,stringsAsFactors = F)#get complete gene list as string
       lgenes=as.vector(lgenes$V1)
       lgenes=tolower(lgenes)
       lgenes=firstup(lgenes)
     }else{lgenes=NULL}
     
     if(input$recgene ==T){
       validate(
         need(input$recgeneli, "Please Upload genelist")
       )
       rgene=input$recgeneli
       rgenes=read.table(rgene$datapath,stringsAsFactors = F)#get complete gene list as string
       rgenes=as.vector(rgenes$V1)
       rgenes=tolower(rgenes)
       rgenes=firstup(rgenes)
     }else{rgenes=NULL}
     
     if(is.null(lgenes)==F & is.null(rgenes)==T ){
       rl=rl[rl$Ligand.ApprovedSymbol %in% lgenes,]
     }else if(is.null(lgenes)==T & is.null(rgenes)==F){
       rl=rl[rl$Receptor.ApprovedSymbol %in% rgenes,]
     }else if(is.null(lgenes)==F & is.null(rgenes)==F ){
       rl=rl[(rl$Receptor.ApprovedSymbol %in% rgenes) & (rl$Ligand.ApprovedSymbol %in% lgenes),]
     }else{rl=rl}
     
     validate(
       need(nrow(rl)>=1, "No Ligand-Receptor Pairs for the combination chosen")
     )
     recid=toupper(rl$Receptor.ApprovedSymbol)
     ligid=toupper(rl$Ligand.ApprovedSymbol)
     urlr= paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",recid,sep = "")
     urll= paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",ligid,sep = "")
     rl$Receptor.ApprovedSymbol=paste0("<a href='",urlr,"'target='_blank'>",rl$Receptor.ApprovedSymbol,"</a>")
     rl$Ligand.ApprovedSymbol=paste0("<a href='",urll,"'target='_blank'>",rl$Ligand.ApprovedSymbol,"</a>")
     return(rl)
   })
   
   
   
   #print data TABLE
   output$ligrecpairs = DT::renderDataTable({
     input$receptor
     input$ligand
     input$recprj
     input$ligprj
     input$rectype
     input$ligtype
     input$liggene
     input$recgene
     input$liggeneli
     input$recgeneli
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       DT::datatable(ligrecpairs(),
                     extensions = c('Buttons','Scroller'),
                     options = list(dom = 'Bfrtip',
                                    searchHighlight = TRUE,
                                    pageLength = 10,
                                    lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                    scrollX = TRUE,
                                    buttons = c('copy', 'print')
                     ),rownames=FALSE,escape = F,selection = list(mode = 'single', selected =1))
     })
   })
   
   output$dwldtab = renderUI({
     downloadButton('dwldtable','Download')
   }) 
   
   output$dwldtable <- downloadHandler(
     filename = function() { paste(input$ligprj,"_",input$ligtype,"_",input$recprj,"_",input$rectype,'.csv', sep='') },
     content = function(file) {
       write.csv(ligrecpairs(), file,row.names=FALSE)
     })
   
   ###################################################
   ###################################################
   ################# kegg pathways ##################
   ###################################################
   ###################################################
   output$rec = DT::renderDataTable({
     input$receptor
     input$ligand
     input$recprj
     input$ligprj
     input$rectype
     input$ligtype
     input$liggene
     input$recgene
     input$liggeneli
     input$recgeneli
     DT::datatable(ligrecpairs(),
                   extensions = c('Buttons','Scroller'),
                   options = list(dom = 'Bfrtip',
                                  searchHighlight = TRUE,
                                  pageLength = 10,
                                  lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                  scrollX = TRUE,
                                  buttons = c('copy', 'print')
                   ),rownames=FALSE,escape = F,selection = list(mode = 'single', selected =1))
   })
   
   keggids <- reactive({
     #get receptor list and annotate it to get entrez ids
     tab=ligrecpairs()
     tab$pair=tab$Pair.Name
     tab= tab %>% separate(pair,c("lig","rec"),sep="_")
     res <- AnnotationDbi::select(org.Mm.eg.db, keys=as.character(tab$rec), columns=c("ENTREZID","SYMBOL"), keytype="SYMBOL")
     res2 <- AnnotationDbi::select(org.Mm.eg.db, keys=as.character(tab$lig), columns=c("ENTREZID","SYMBOL"), keytype="SYMBOL")
     res <- subset(res,!duplicated(res$ENTREZID))
     res2 <- subset(res2,!duplicated(res2$ENTREZID))
     final_tab<- left_join(tab,res, by=c("rec"="SYMBOL")) %>% rename("ENTREZID"="ENTREZID_rec")
     final_tab<- left_join(final_tab,res2, by=c("lig"="SYMBOL")) %>% rename("ENTREZID"="ENTREZID_lig")
     final_tab=unique(final_tab)
     # final_tab$ENTREZID_rec=paste0("ncbi-geneid:",final_tab$ENTREZID_rec,sep="")
     # final_tab$ENTREZID_lig=paste0("ncbi-geneid:",final_tab$ENTREZID_lig,sep="")
     return(final_tab)
   })
   gageres <- reactive({
     final_tab=keggids()
     s = input$rec_rows_selected #select rows from table
     select_tab = final_tab[s, , drop=FALSE]
     select_tab$ENTREZID_rec=paste0("ncbi-geneid:",select_tab$ENTREZID_rec,sep="")
     keggids=(keggConv("mmu",select_tab$ENTREZID_rec))
     
     #Find pathways
     res=keggLink("pathway",keggids)
     res2=as.data.frame(res)
     res2$Receptor_id=names(res)
     # res3=res2[res2$res %in% unique(res2$res[duplicated(res2$res)]),]
     # table=as.data.frame(table(res3$res))
     # table=table[order(-table$Freq),]
     table=res2 %>% tidyr::separate(res,c("path","pathway_id")) %>% dplyr::select(-path)
     # final_tab=final_tab %>% separate(ENTREZID,c("ncbi","ENTREZID"),sep=":") %>% dplyr::select(-ncbi)
     for(i in 1: nrow(table)){
       table$Name[i]=keggGet(table$pathway_id[i])[[1]]$PATHWAY_MAP
       # genes=keggGet(table$pathway_id[i])[[1]]$GENE
       # ind=seq(1,length(genes),2)
       # genes_entrez=genes[ind]
       # genelist=final_res$ENTREZID[final_res$ENTREZID %in% genes_entrez]
       allgenelist=keggLink("mmu",table$pathway_id[i]) #for each kegg id, get gene list
       p=strsplit(allgenelist,":")
       genes_entrez=sapply(p,"[",2)
       genelist=genes_entrez[genes_entrez %in% final_tab$ENTREZID_rec]
       genelist=unique(genelist)
       table$Num_of_Rec_genes_in_pathway[i]=length(genelist)
     }
     return(table)
   })
   
   output$Keggpaths = DT::renderDataTable({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       DT::datatable(gageres(),
                     extensions = c('Buttons','Scroller'),
                     options = list(dom = 'Bfrtip',
                                    searchHighlight = TRUE,
                                    pageLength = 10,
                                    lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                    scrollX = TRUE,
                                    buttons = c('copy', 'print')
                     ),rownames=FALSE,escape = F,selection = list(mode = 'single', selected =1))
     })
   })
   ###################################################
   ###################################################
   ################ PATHWAY ANALYSIS #################
   ###################################################
   ###################################################
   output$plots = renderImage({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       #input$makeplot
       path=gageres() #get KEGG id's  of pathways
       s = input$Keggpaths_rows_selected #select rows from table
       path = path[s, , drop=FALSE]#get data corresponding to selected row in table
       pId=path$pathway_id
       path=path %>% tidyr::separate(Receptor_id,c("rec","rec_id")) %>% dplyr::select(-rec)
       recid=unique(path$rec_id)
       
       final_res=keggids()
       allgenelist=keggLink("mmu",pId) #for each kegg id, get gene list
       p=strsplit(allgenelist,":")
       genes_entrez=sapply(p,"[",2)
       rec_genes=final_res$ENTREZID_rec
       final_res=final_res[final_res$ENTREZID_rec %in% recid,]
       allgenes=unique(c(rec_genes,final_res$ENTREZID_lig))
       genelist=genes_entrez[genes_entrez %in% allgenes]
       genelist=unique(genelist)
       genelist=paste0("mmu:",genelist,sep="")
       
       myurl=mark.pathway.by.objects(pId,genelist) #get url of pathway image
       outfile = tempfile(fileext='.png') #create temp file
       png(outfile, width=900, height=800) #get temp file in png format
       download.file(myurl,outfile,mode="wb") #download png into the temp file
       png = readPNG(outfile) # read the PNG from the temp file and display
       dev.off()
       
       list(src = outfile,contentType = 'image/png',width = 900,height = 800,alt = "This is alternate text")
     })
   }, deleteFile = TRUE) #delete temp file
   
}#end of server