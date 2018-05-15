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
library(data.table)
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
      plot1=TSNEPlot(object = scrna,group.by = "ident",no.legend = FALSE,do.label = TRUE, do.return=T, pt.size = input$pointa) + theme(legend.position="bottom")
    }else if(input$categorya =="var" & input$tsnea %in% tsne){
      plot1=TSNEPlot(object = scrna,group.by = tsnea,no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$pointa)
    }else if(input$categorya =="var" & input$tsnea %in% feature){
      plot1=FeaturePlot(object = scrna, features.plot = tsnea, cols.use = c("grey", "blue"),reduction.use = "tsne",do.return=T,pt.size = input$pointa)
      plot1=eval(parse(text=paste("plot1$`",tsnea,"`",sep="")))
    }

    markers=markergenes()
      s=input$markergenes_rows_selected # get  index of selected row from table
      markers=markers[s, ,drop=FALSE]
      plot2=FeaturePlot(object = scrna, features.plot = rownames(markers), cols.use = c("grey","blue"),reduction.use = "tsne",
                        no.legend = FALSE,pt.size = input$pointa,do.return = T)
      plot2=eval(parse(text=paste("plot2$",rownames(markers),sep="")))
    
      plot3=VlnPlot(object = scrna, features.plot = rownames(markers),group.by = input$grptype,do.return = T,x.lab.rot=TRUE)
      plot4=RidgePlot(object = scrna, features.plot = rownames(markers),group.by = input$grptype,do.return = T,x.lab.rot=TRUE)
      
    
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
  
  output$identdef = renderUI({
    scrna=fileload()
    options=names(scrna@misc)
    selectInput("identdef", "First cluster/variable of comparison",options)
  })
  
  output$identa = renderUI({
    scrna=fileload()
    if(input$setident==T){
      scrna <- SetAllIdent(object = scrna, id = input$setidentlist)
      options=unique(scrna@ident)
    }else{
      options=levels(scrna@ident)
    }
    selectInput("identa", "First cluster to compare",options)
    })
  
    output$identb = renderUI({
      scrna=fileload()
      if(input$setident==T){
        scrna <- SetAllIdent(object = scrna, id = input$setidentlist)
        options=unique(scrna@ident)
      }else{
        options=levels(scrna@ident)
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
    }
    if(input$setident==F){
      markers=eval(parse(text=paste("scrna@misc$`",input$identdef,"`",sep="")))
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
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
     scrna=fileload()
     options=levels(scrna@ident)
     selectInput("clust1","Pick cluster1",options,selected=options[1])
     })
   })
   
   output$clust2 <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
     scrna=fileload()
     options=levels(scrna@ident)
     selectInput("clust2","Pick cluster2",options, selected=options[2])
     })
   })
   
   output$list1.1 <- renderUI({
     fileInput('genelist1.1', 'Upload Gene List1',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
   })
   
   output$list2.1 <- renderUI({
     fileInput('genelist2.1', 'Upload Gene List2',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
   })
   
   output$clust1.1 <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
     scrna=fileload()
     options=levels(scrna@ident)
     selectInput("clust1.1","Pick cluster1",options,selected=options[1])
     })
   })
   
   output$clust2.1 <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
     scrna=fileload()
     options=levels(scrna@ident)
     selectInput("clust2.1","Pick cluster2",options,selected=options[2])
     })
   })
   ###################################################
   ###################################################
   #### Load lig-receptor list and display results ###
   ###################################################
   ###################################################
   
   datasetInput = reactive({
     scrna=fileload()
     tt=rownames(scrna@raw.data)
     file = read.csv("data/param.csv")
     org=as.character(file$organism[file$projects==input$projects])
     genes=fread("data/ligrecgenes.txt",header = TRUE)
     if(org=="human"){
       genes$genes=toupper(genes$genes)
     }
     genes2=tt[tt %in% genes$genes]
     my.data=FetchData(scrna,c("ident","nGene",genes2))
     my.data= my.data %>% rename('clust'='ident')
     
     
     rl=read.csv("data/lig-rec.csv")
     if(org=="human"){
       rl= rl %>% dplyr::select(Pair.Name:Receptor.ApprovedSymbol)
       rl$Pair.Name=toupper(rl$Pair.Name)
       rl$Ligand.ApprovedSymbol=toupper(rl$Ligand.ApprovedSymbol)
       rl$Receptor.ApprovedSymbol=toupper(rl$Receptor.ApprovedSymbol)
     }else if(org=="mouse"){
       rl= rl %>% dplyr::select(Mouse_LigandSym:Mouse.Pairs) %>% rename("Pair.Name"="Mouse.Pairs","Ligand.ApprovedSymbol"="Mouse_LigandSym","Receptor.ApprovedSymbol"="Mouse_RecSym")
     }
     result=data.frame()
     res=data.frame()
     for(i in 1:(length(unique(my.data$clust)))){
       for(j in 1:(length(unique(my.data$clust)))){
         if(i!=j){
           test=my.data[my.data$clust==levels(my.data$clust)[i] | my.data$clust==levels(my.data$clust)[j],]
           R_c1=test[test$clust==levels(my.data$clust)[i] ,(colnames(test) %in% rl$Receptor.ApprovedSymbol)]
           L_c2=test[test$clust==levels(my.data$clust)[j] , (colnames(test) %in% rl$Ligand.ApprovedSymbol)]
           keep1 = colSums(R_c1>1)>=.5*dim(R_c1)[1]
           keep2 = colSums(L_c2>1)>=.5*dim(L_c2)[1]
           
           R_c1=R_c1[,keep1]
           L_c2=L_c2[,keep2]
           res=rl[(rl$Ligand.ApprovedSymbol %in% colnames(L_c2)) & (rl$Receptor.ApprovedSymbol %in% colnames(R_c1)),]
           
         }
         else{}
         if(nrow(res)!=0){
           res$Receptor_cluster=levels(my.data$clust)[i]
           res$Lig_cluster=levels(my.data$clust)[j]
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
   
}#end of server