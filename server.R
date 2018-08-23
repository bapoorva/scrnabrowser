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
cpallette=c("#64B2CE", "#DA5724", "#74D944", "#CE50CA", "#C0717C", "#CBD588", "#5F7FC7",
            "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
            "#8A7C64", "#599861","#000099","#FFCC66","#99CC33","#CC99CC","#666666")

my_username <- c("Sealelab","Morriseylab","Jainlab","allusers")
my_password <- c("pseale@999#","emorrisey$123","rjain@2018","allusers@1")

server <- function(input, output,session) {
  
  ###################################################
  ###################################################
  ####### Display project list and load data  #######
  ###################################################
  ###################################################
   values <- reactiveValues(authenticated = FALSE)
  
   # Return the UI for a modal dialog with data selection input. If 'failed'
   # is TRUE, then display a message that the previous value was invalid.
   dataModal <- function(failed = FALSE) {
     modalDialog(
       textInput("username", "Username:"),
       passwordInput("password", "Password:"),
       footer = tagList(
         # modalButton("Cancel"),
         actionButton("ok", "OK")
       )
     )
   }
   
   # Show modal when button is clicked.
   # This `observe` is suspended only whith right user credential
   
   obs1 <- observe({
     showModal(dataModal())
   })
  
  # When OK button is pressed, attempt to authenticate. If successful,
  # remove the modal.
  obs2 <- observe({
    req(input$ok)
    isolate({
      Username <- input$username
      Password <- input$password
    })
    Id.username <- which(my_username == Username)
    Id.password <- which(my_password == Password)
    if (length(Id.username) > 0 & length(Id.password) > 0) {
      if (Id.username == Id.password) {
        Logged <<- TRUE
        values$authenticated <- TRUE
        obs1$suspend()
        removeModal()

      } else {
        values$authenticated <- FALSE
      }
    }
  })
  # output$dataInfo <- renderPrint({
  #   if (values$authenticated) "OK!!!!!"
  #   else "You are NOT authenticated"
  # })
  
  #Read the parameter file
  readexcel = reactive({
     user=input$username
     file = read.csv(paste("data/",user,"_param.csv",sep=""))
    #file = read.csv("data/Morriseylab_param.csv")
  })
  
  #Get Project list and populate drop-down
  output$projects = renderUI({
    excel=readexcel()
    prj=as.character(excel$projects)
    selectInput("projects","Select a project",as.list(sort(as.character(prj))))
  })
  
  #display project list in Dashboard
  output$datasetTable<- renderTable({
    user=input$username
    read.csv(paste("data/",user,"_param.csv",sep=""))
    #read.csv('data/param.csv',stringsAsFactors = F)
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

    if(input$categorya2 =="clust" & input$subsa==F){
      plot1=DimPlot(object = scrna,reduction.use=input$umapa,group.by = "ident",no.legend = FALSE,do.label = TRUE, do.return=T, pt.size = input$pointa2,label.size = 7, cols.use=cpallette)
    }else if(input$categorya2 =="clust" & input$subsa==TRUE){
      cells=names(scrna@ident[scrna@ident==input$selclust])
      plot1=DimPlot(object = scrna,reduction.use=input$umapa,cells.highlight=cells,group.by = "ident",no.legend = FALSE,do.label = F, do.return=T, pt.size = input$pointa2, cols.use=cpallette)
    }else if(input$categorya2=="geneexp"){
      plot1=FeaturePlot(object = scrna,reduction.use=input$umapa, features.plot = input$gene1a, cols.use = c("grey", "blue"),do.return=T,pt.size = input$pointa2,no.legend = FALSE)
      plot1=eval(parse(text=paste("plot1$",input$gene1a,sep="")))
    }else if(input$categorya2 =="var" & input$tsnea2 %in% tsne & input$subsa==FALSE){
      plot1=DimPlot(object = scrna,reduction.use=input$umapa,group.by = tsnea,no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$pointa2,label.size = 7, cols.use=cpallette)
    }else if(input$categorya2 =="var" & input$tsnea2 %in% tsne & input$subsa==TRUE){
      t=paste("rownames(scrna@meta.data[scrna@meta.data$",input$tsnea2,"==\"",input$selclust2,"\",])",sep="")
      cells=eval(parse(text=t))
      plot1=DimPlot(object = scrna,reduction.use=input$umapa,group.by = tsnea,cells.highlight=cells,no.legend = FALSE,do.label =F, do.return=T,pt.size = input$pointa2, cols.use=cpallette)
    }else if(input$categorya2 =="var" & input$tsnea2 %in% feature & input$subsa==FALSE){
      plot1=FeaturePlot(object = scrna,reduction.use=input$umapa, features.plot = tsnea, cols.use = c("grey", "blue"),do.return=T,pt.size = input$pointa2,no.legend = FALSE)
      plot1=eval(parse(text=paste("plot1$",tsnea,sep="")))
    }else if(input$categorya2 =="var" & input$tsnea2 %in% feature & input$subsa==TRUE){
      t=paste('rownames(scrna@meta.data[scrna@meta.data$',input$tsnea2, '>',input$tsnea2lim[1], ' & metadata$',input$tsnea2, '<', input$tsnea2lim[2],',])',sep="")
      cells=eval(parse(text=t))
      plot1=FeaturePlot(object = scrna,reduction.use=input$umapa, features.plot = tsnea,cells.use = cells, cols.use = c("grey", "blue"),do.return=T,pt.size = input$pointa2,no.legend = FALSE)
      plot1=eval(parse(text=paste("plot1$",tsnea,sep="")))
     }

    if(input$categoryb2 =="clust" & input$subsb==F){
      plot2=DimPlot(object = scrna,reduction.use=input$umapb,group.by = "ident",no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$pointa2,label.size = 7, cols.use=cpallette)
    }else if(input$categoryb2 =="clust" & input$subsb==TRUE){
      cells=names(scrna@ident[scrna@ident==input$selclustb])
      plot2=DimPlot(object = scrna,reduction.use=input$umapb,cells.highlight=cells,group.by = "ident",no.legend = FALSE,do.label = F, do.return=T, pt.size = input$pointa2, cols.use=cpallette)
    }else if(input$categoryb2=="geneexp"){
      plot2=FeaturePlot(object = scrna,reduction.use=input$umapb, features.plot = input$gene2a, cols.use = c("grey", "blue"),do.return=T,pt.size = input$pointa2,no.legend = FALSE)
      plot2=eval(parse(text=paste("plot2$",input$gene2a,sep="")))
    }else if(input$categoryb2 =="var" & input$tsneb2 %in% tsne & input$subsb==F){
      plot2=DimPlot(object = scrna,reduction.use=input$umapb,group.by = tsneb,no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$pointa2,label.size = 7, cols.use=cpallette)
    }else if(input$categoryb2 =="var" & input$tsneb2 %in% tsne & input$subsb==TRUE){
      t=paste("rownames(scrna@meta.data[scrna@meta.data$",input$tsneb2,"==\"",input$selclustb2,"\",])",sep="")
      cells=eval(parse(text=t))
      plot2=DimPlot(object = scrna,reduction.use=input$umapb,group.by = tsneb,cells.highlight=cells,no.legend = FALSE,do.label = F, do.return=T,pt.size = input$pointa2, cols.use=cpallette)
    }else if(input$categoryb2 =="var" & input$tsneb2 %in% feature & input$subsb==F){
      plot2=FeaturePlot(object = scrna,reduction.use=input$umapb, features.plot = tsneb, cols.use = c("grey", "blue"),do.return=T,pt.size = input$pointa2,no.legend = FALSE)
      plot2=eval(parse(text=paste("plot2$",tsneb,sep="")))
    }else if(input$categoryb2 =="var" & input$tsneb2 %in% feature & input$subsb==TRUE){
      t=paste('rownames(scrna@meta.data[scrna@meta.data$',input$tsneb2, '>',input$tsneb2lim[1], ' & metadata$',input$tsneb2, '<', input$tsneb2lim[2],',])',sep="")
      cells=eval(parse(text=t))
      plot2=FeaturePlot(object = scrna,reduction.use=input$umapb, features.plot = tsneb,cells.use = cells, cols.use = c("grey", "blue"),do.return=T,pt.size = input$pointa2,no.legend = FALSE)
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
      paste0("Compare_tsne.pdf")
    },
    content = function(file){
      pdf(file,width=14,height = 8,useDingbats=FALSE)
      plot(comptsne2())
      dev.off()
    })
  
  ###################################################
  ###################################################
  ########### Interactive Tsne plots  ###############
  ###################################################
  ###################################################
  output$intervar = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data) 
    #metadata=metadata %>% select(starts_with("var"))
    var=colnames(metadata)
    selectInput("intervar","Select a Variable",var,"pick one")
  })
  
  output$setcategory = renderUI({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    metadata=metadata %>% select(starts_with("var_"))
    var=c(colnames(metadata))
    selectInput("setcategory","Choose category",var,"pick one")
  })
  
  intertsne = reactive({
    scrna=fileload()
    plot1=DimPlot(object = scrna,reduction.use=input$umapint,group.by = input$setcategory,no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$umap_pointsize,label.size = 5, cols.use=cpallette)
    plot=ggplotly(plot1)
    return(plot)
  })
  
  output$intertsne = renderPlotly({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      intertsne()
    })
  })
  
  intergene = reactive({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    met= sapply(metadata,is.numeric)
    #metadata=metadata %>% select(starts_with("var"))
    tsnea=input$intervar
    feature=names(met[met==TRUE])
    #feature=c("nGene","nUMI","percent.mito","S.Score","G2M.Score","var.ratio.pca")
    tsne=names(met[met==FALSE])
    
    if(input$intercat=="geneexp"){
      plot1=FeaturePlot(object = scrna,reduction.use=input$umapint, features.plot = input$geneinter, cols.use = c("grey", "blue"),do.return=T,pt.size = input$umap_pointsize,no.legend = FALSE)
      plot1=eval(parse(text=paste("plot1$",input$geneinter,sep="")))
    }else if(input$intercat =="var" & tsnea %in% tsne){
      plot1=DimPlot(object = scrna,reduction.use=input$umapint,group.by = tsnea,no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$umap_pointsize,label.size = 7, cols.use=cpallette)
    }else if(input$intercat =="var" & tsnea %in% feature){
      plot1=FeaturePlot(object = scrna,reduction.use=input$umapint, features.plot = tsnea, cols.use = c("grey", "blue"),do.return=T,pt.size = input$umap_pointsize,no.legend = FALSE)
      plot1=eval(parse(text=paste("plot1$",tsnea,sep="")))
    }
    plot=ggplotly(plot1)
    return(plot)
  })
  
  output$intergene = renderPlotly({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      intergene()
    })
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
    var=c(colnames(metadata),'Cell.group')
    selectInput("tsnea","Select tSNE group to display",var,selected = "Cell.group")
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
    scrna@meta.data$var_cluster=as.numeric(as.character(scrna@meta.data$var_cluster))
    tsnea=input$tsnea
    tsneb=input$tsneb
    feature=names(met[met==TRUE])
    tsne=names(met[met==FALSE])
    
    if(input$tsnea =="Cell.group"){
      plot1=DimPlot(object = scrna,reduction.use=input$umapdeg,group.by = "ident",no.legend = FALSE,do.label = TRUE, do.return=T, pt.size = input$pointa,label.size = 7,cols.use=cpallette) + theme(legend.position="bottom")
    }else if(input$tsnea %in% tsne){
      plot1=DimPlot(object = scrna,reduction.use=input$umapdeg,group.by = tsnea,no.legend = FALSE,do.label = TRUE, do.return=T,pt.size = input$pointa,label.size = 7,cols.use=cpallette) + theme(legend.position="bottom")
    }else if(input$tsnea %in% feature){
      plot1=FeaturePlot(object = scrna, features.plot = tsnea, cols.use = c("grey", "blue"),reduction.use = input$umapdeg,do.return=T,pt.size = input$pointa)
      plot1=eval(parse(text=paste("plot1$`",tsnea,"`",sep="")))
    }
   
    markers=markergenes()
      s=input$markergenes_rows_selected # get  index of selected row from table
      markers=markers[s, ,drop=FALSE]
      plot2=FeaturePlot(object = scrna, features.plot = rownames(markers), cols.use = c("grey","blue"),reduction.use = input$umapdeg,
                        no.legend = FALSE,pt.size = input$pointa,do.return = T)
      plot2=eval(parse(text=paste("plot2$`",rownames(markers),"`",sep="")))
      if(input$checkviolin ==T){
      plot3=VlnPlot(object = scrna, features.plot = rownames(markers),group.by = input$setidentlist,do.return = T,x.lab.rot=TRUE,point.size.use=NA,cols.use=cpallette)
      }else{plot3=VlnPlot(object = scrna, features.plot = rownames(markers),group.by = input$setidentlist,do.return = T,x.lab.rot=TRUE,cols.use=cpallette)}
      plot4=RidgePlot(object = scrna, features.plot = rownames(markers),group.by = input$setidentlist,do.return = T,x.lab.rot=TRUE,cols.use=cpallette)
      
    
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
      paste0("Plot.pdf")
    },
    content = function(file){
      pdf(file, width = 12, height = 11,useDingbats=FALSE)
      plot(comptsne())
      dev.off()
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
      #textInput("bigene_genea", label = "Gene A",value = bigene_genea)
    r<-getGeneRange(fileload(),input$bigene_genea)
    sliderInput("bigene_rangea", "Expression Limit Gene A (log2 UMI)",
                min = 0, max = r[2], value = c(r[1],r[2]),step=.25)
  })
  })
  
  output$bigene_rangeb <- renderUI({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      #textInput("bigene_geneb", label = "Gene B",value = bigene_geneb)
    r<-getGeneRange(fileload(),input$bigene_geneb)
    sliderInput("bigene_rangeb", "Expression Limit Gene B (log2 UMI)",
                min = 0, max = r[2], value = c(r[1],r[2]),step=.25)
  })
  })
  
  output$bigene_rangea2 <- renderUI({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      table=finalres()
      s=input$pairs_res_rows_selected
      table=table[s, ,drop=FALSE]
      validate(need(nrow(table) > 0,"No Ligand-receptor pairs found"))
      bigene_genea=table$ligand
      #textInput("bigene_genea", label = "Gene A",value = bigene_genea)
      r<-getGeneRange(fileload(),bigene_genea)
      sliderInput("bigene_rangea", "Expression Limit of Ligand Gene (log2 UMI)",
                  min = 0, max = r[2], value = c(r[1],r[2]),step=.25)
    })
  })
  
  output$bigene_rangeb2 <- renderUI({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      table=finalres()
      validate(need(nrow(table) > 0,"No Ligand-receptor pairs found"))
      s=input$pairs_res_rows_selected
      table=table[s, ,drop=FALSE]
      bigene_geneb=table$receptor
      #textInput("bigene_geneb", label = "Gene B",value = bigene_geneb)
      r<-getGeneRange(fileload(),bigene_geneb)
      sliderInput("bigene_rangeb", "Expression Limit of Receptor Gene (log2 UMI)",
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

  output$bigeneplot2 <- renderPlot({
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      table=finalres()
      s=input$pairs_res_rows_selected
      table=table[s, ,drop=FALSE]
      validate(need(nrow(table) > 0,"No Ligand-receptor pairs found"))
      bigene_genea=as.character(table$ligand)
      bigene_geneb=as.character(table$receptor)
      bigene_plot(fileload(),
                  c(bigene_genea,bigene_geneb),
                  limita=input$bigene_rangea,
                  limitb=input$bigene_rangeb,
                  marker_size = input$bigene_pointsize2)
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
    selectInput("identa", "First Cell group to compare",options)
    })
  
    output$identb = renderUI({
      scrna=fileload()
      
      if(input$setident==T){
        scrna <- SetAllIdent(object = scrna, id = input$setidentlist)
        options=unique(scrna@ident)
      }else{
        options=levels(scrna@ident)
      }
      checkboxGroupInput("identb", label="Second Cell group to compare",choices=options)
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
      if(input$goButton == 0)
        return()
      isolate({
      scrna <- SetAllIdent(object = scrna, id = input$setidentlist)
    validate(need(input$identb,"Select at least one option from Second cell group to compare to first cell group. If you want to compare to all, uncheck the 'Check to choose a different category to compare' option"))
    validate(need(input$identb!=input$identa,"First and second cell groups can't be the same"))
    
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
    }
    if(input$setident==F){
      markers=eval(parse(text=paste("scrna@misc$`",input$identdef,"`",sep="")))
      geneid=rownames(markers)
      url= paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",geneid,sep = "")
      markers$Link=paste0("<a href='",url,"'target='_blank'>",rownames(markers),"</a>")
    }
    
    })
    return(markers)
  })
  
  output$markergenes = DT::renderDataTable({
    #input$identa
    #input$identb
    input$goButton
    DT::datatable(markergenes(),
                  extensions = 'Buttons', options = list(
                    dom = 'Bfrtip',
                    pageLength = 10,
                    lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                    buttons = list()),
                  rownames=TRUE,selection = list(mode = 'single', selected =1),escape=FALSE)
   
  })
  
  output$downloaddeg <- downloadHandler(
    filename = function() { paste(input$projects, '.csv', sep='') },
    content = function(file) {
      write.csv(markergenes(), file)
    })
  

  
  # ####################################################
  # ###################################################
  # ####### Display Violin  plot with controls  #######
  # ###################################################
  # ###################################################
  # output$grptype = renderUI({
  #   withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
  #   scrna=fileload()
  #   metadata=as.data.frame(scrna@meta.data)
  #   metadata=metadata %>% select(starts_with("var_"))
  #   var=c("ident",colnames(metadata))
  #   selectInput("grptype","Select a Variable to group the Violin Plot",var,"pick one")
  #   })
  # })
  
  
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
       metadata=metadata %>% select(starts_with("var_"))
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
       paste0("Heatmap.pdf")
     },
     content = function(file){
       pdf(file, width = 13, height = 8,useDingbats=FALSE)
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
   
   output$pairby <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
       scrna=fileload()
       metadata=as.data.frame(scrna@meta.data)
       metadata=metadata %>% select(starts_with("var_"))
       options=colnames(metadata)
       selectInput("pairby","Select cell group ",options,selected=options[1])
     })
   })
   
   output$clust1 <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
     scrna=fileload()
     t=paste("scrna@meta.data$",input$pairby,sep="")
     options=unique(eval(parse(text=t)))
     if(input$pairby=="ident"){options=levels(scrna@ident)}
     selectInput("clust1","Pick cellgroup for Receptor",options,selected=options[1])
     })
   })
   
   output$clust2 <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
     scrna=fileload()
     t=paste("scrna@meta.data$",input$pairby,sep="")
     options=unique(eval(parse(text=t)))
     if(input$pairby=="ident"){options=levels(scrna@ident)}
     selectInput("clust2","Pick cellgroup for Ligand",options, selected=options[2])
     })
   })
   
   output$list1.1 <- renderUI({
     fileInput('genelist1.1', 'Upload Receptor Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
   })
   
   output$list2.1 <- renderUI({
     fileInput('genelist2.1', 'Upload Ligand Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
   })
   
   output$clust1.1 <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
     scrna=fileload()
     t=paste("scrna@meta.data$",input$pairby2,sep="")
     options=unique(eval(parse(text=t)))
     selectInput("clust1.1","Pick cellgroup for Receptor",options,selected=options[1])
     })
   })
   
   output$clust2.1 <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
     scrna=fileload()
     t=paste("scrna@meta.data$",input$pairby2,sep="")
     options=unique(eval(parse(text=t)))
     selectInput("clust2.1","Pick cellgroup for Ligand",options,selected=options[2])
     })
   })
   ###################################################
   ###################################################
   #### Load lig-receptor list and display results ###
   ###################################################
   ###################################################
   output$source <- renderUI({
     file = read.csv("data/param.csv")
     org=as.character(file$organism[file$projects==input$projects])
     if(org=="mouse"){rl=read.csv("data/Mm_PairsLigRec.csv")}else if(org=="human"){rl=read.csv("data/Hs_PairsLigRec.csv")}
     options=as.character(unique(rl$Pair.Source))
     checkboxGroupInput('source', label='Select source(s)',choices=options,selected=options[1])
   })
   
   output$evidence <- renderUI({
     file = read.csv("data/param.csv")
     org=as.character(file$organism[file$projects==input$projects])
     if(org=="mouse"){rl=read.csv("data/Mm_PairsLigRec.csv")}else if(org=="human"){rl=read.csv("data/Hs_PairsLigRec.csv")}
     options=as.character(unique(rl$Pair.Evidence))
     checkboxGroupInput('evidence',label='Select Evidence(s)',choices=options,selected=options[1])
   })
   
   datasetInput = reactive({
     scrna=fileload()
     var=as.character(input$pairby)
     tt=rownames(scrna@raw.data)
     file = read.csv("data/param.csv")
     org=as.character(file$organism[file$projects==input$projects])
     genes=fread("data/ligrecgenes.txt",header = TRUE)
     if(org=="human"){
       genes$genes=toupper(genes$genes)
     }
     genes2=tt[tt %in% genes$genes]
     my.data=FetchData(scrna,c(var,"nGene",genes2))
     colnames(my.data)[1]= "clust"
     #my.data$clust=factor(my.data$clust,levels=unique(my.data$clust))
     
     if(org=="mouse"){rl=read.csv("data/Mm_PairsLigRec.csv")}else if(org=="human"){rl=read.csv("data/Hs_PairsLigRec.csv")}
     result=data.frame()
     res=data.frame()
     for(i in 1:(length(levels(my.data$clust)))){
       for(j in 1:(length(levels(my.data$clust)))){
          if(i!=j){
           test=my.data[my.data$clust==levels(my.data$clust)[i] | my.data$clust==levels(my.data$clust)[j],]
           R_c1=test[test$clust==levels(my.data$clust)[i] ,(colnames(test) %in% rl$receptor)]
           L_c2=test[test$clust==levels(my.data$clust)[j] , (colnames(test) %in% rl$ligand)]
           keep1 = colSums(R_c1>1)>=.5*dim(R_c1)[1]
           keep2 = colSums(L_c2>1)>=.5*dim(L_c2)[1]
           R_c1=R_c1[,keep1]
           L_c2=L_c2[,keep2]
           res=rl[(rl$ligand %in% colnames(L_c2)) & (rl$receptor %in% colnames(R_c1)),]

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
     return(result)
   })
   
   finalres= reactive({
     result=datasetInput()
     if(input$clust=="clust" & input$gene=="allgene"){
       #clusters=c(input$clust1,input$clust2)
       result=result[(result$Receptor_cluster %in% input$clust1) & (result$Lig_cluster%in% input$clust2),]
     }else if(input$clust=="clust" & input$gene=="genelist"){
       #clusters.1=c(input$clust1.1,input$clust2.1)
       result=result[(result$Receptor_cluster %in% input$clust1.1) & (result$Lig_cluster%in% input$clust2.1),]
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
       result=result[(result$receptor %in% g1) & (result$ligand %in% g2),]
     }else{
       result=result
     }
     if(input$checksource==T){result=result[result$Pair.Source %in% input$source,]}
     if(input$checkevi==T){result=result[result$Pair.Evidence %in% input$evidence,]}
     return(result)
   })
   
   #print data TABLE
   output$pairs_res = DT::renderDataTable({
     input$pairby
     input$clust1
     input$clust2
     input$genelist1
     input$genelist2
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
       DT::datatable(finalres(),
                     extensions = c('Buttons','Scroller'),
                     options = list(dom = 'Bfrtip',
                                    searchHighlight = TRUE,
                                    pageLength = 10,
                                    lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                    scrollX = TRUE,
                                    buttons = c('copy', 'print')
                     ),rownames=FALSE,caption= "Result",selection = list(mode = 'single', selected =1),escape = F)
     })
   })
   
   ###################################################
   ###################################################
   ########### Plot gene expression  ################
   ###################################################
   ###################################################
   
   geplots = reactive({
     scrna=fileload()
     validate(need(input$geneid,"Enter the gene symbol"))
     plot2=FeaturePlot(object = scrna, features.plot = input$geneid, cols.use = c("grey","blue"),reduction.use = input$umapge,
                       no.legend = FALSE,pt.size = input$genenid_pointsize,do.return = T)
     plot2=eval(parse(text=paste("plot2$`",input$geneid,"`",sep="")))
     plot3=VlnPlot(object = scrna, features.plot = input$geneid,group.by = "ident",do.return = T,x.lab.rot=TRUE,cols.use=cpallette)
     plot4=RidgePlot(object = scrna, features.plot = input$geneid,group.by = "ident",do.return = T,x.lab.rot=TRUE,cols.use=cpallette)
     
     
     row1=plot_grid(plot2,align = 'h', rel_heights = c(1, 1),axis="lr", nrow=1)
     row2=plot_grid(plot3,plot4,align = 'h', rel_heights = c(1, 1),axis="lr", nrow=1)
     plot_grid(row1,row2,align = 'v', rel_heights = c(1.7, 1),axis="tb",ncol=1)
     
   })
   
   output$geplots = renderPlot({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       geplots()
     })
   })
   
   output$downloadplotge <- downloadHandler(
     filename = function() {
       paste0(input$geneid,"_Geneexp_plot.pdf",sep="")
     },
     content = function(file){
       pdf(file,width=8,height = 13,useDingbats=FALSE)
       plot(geplots())
       dev.off()
     })
   
   ###################################################
   ###################################################
   ########### Plot gene expression in clusters  #####
   ###################################################
   ###################################################
   output$setvar = renderUI({
     scrna=fileload()
     metadata=as.data.frame(scrna@meta.data)
     metadata=metadata %>% select(starts_with("var_"))
     var=c(colnames(metadata))
     selectInput("setvar","Choose category",var,"pick one")
   })
   
   output$selectcluster = renderUI({
     scrna=fileload()
       options <- paste0("scrna@meta.data$",input$setvar,sep="")
       options=unique(eval(parse(text=options)))
       options=options[order(options)]
     selectInput("selectcluster", "Select Cell group",options)
   })
   
   clusts= reactive({
     scrna=fileload()
     scrna <- SetAllIdent(object = scrna, id = input$setvar)
     avgexp=AverageExpression(object = scrna)
     avgexp= avgexp %>% dplyr::select(input$selectcluster)
     genes.use=rownames(avgexp)
     data.use <- GetAssayData(object = scrna,slot = "data")
     cells <- WhichCells(object = scrna, ident = input$selectcluster)
     thresh.min=0
     data.temp <- round(x = apply(X = data.use[genes.use, cells, drop = F],
         MARGIN = 1,
         FUN = function(x) {
           return(sum(x > thresh.min) / length(x = cells))
         }),digits = 3)
     names(data.temp)=genes.use
     data.temp=as.data.frame(data.temp)
     data.temp$id=rownames(data.temp)
     avgexp$id=rownames(avgexp)
     df=inner_join(avgexp,data.temp,by="id") 
     rownames(df)=df$id
     df= df %>% dplyr::select(input$selectcluster,data.temp) 
     df=df[order(-df[,1],-df[,2]),]
     colnames(df)= c("Average Expression","Percentage of cells expressed in")
     df$max_avg=max(df$`Average Expression`)
     df$min_avg=min(df$`Average Expression`)
     return(df)
   })
   
   clustable= reactive({
     df=clusts()
     df= df %>% dplyr::select(-max_avg:-min_avg)
     df=df[df$`Percentage of cells expressed in` >input$pctslider[1] & df$`Percentage of cells expressed in` <input$pctslider[2],]
     df=df[df$`Average Expression` >input$avgexpslider[1] & df$`Average Expression` <input$avgexpslider[2],]
     return(df)
   })
   
   output$pctslider <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
       df=clusts()
       sliderInput("pctslider", "Percent Expressed:",min =0, max = 1, value =c(0.1,1))
     })
   })
   
   output$avgexpslider <- renderUI({
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
       df=clusts()
       min=unique(df$min_avg)
       max=unique(df$max_avg)
       sliderInput("avgexpslider", "Average Expression:",min = min, max = max, value = c(min,max))
     })
   })
   
   output$clustable = DT::renderDataTable({
     input$setvar
     input$selectcluster
     input$avgexpslider
     input$pctslider
     input$umapclust
     withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
       DT::datatable(clustable(),
                     extensions = c('Buttons','Scroller'),
                     options = list(dom = 'Bfrtip',
                                    searchHighlight = TRUE,
                                    pageLength = 10,
                                    lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                    scrollX = TRUE,
                                    buttons = c('copy', 'print')
                     ),rownames=TRUE,caption= "Cluster-wise Gene expression",selection = list(mode = 'single', selected =1),escape = F)
     })
   })
   
   
   clustplots= reactive({
     scrna=fileload()
     tab=clustable()
     s=input$clustable_rows_selected
     tab=tab[s, ,drop=FALSE]
     gene=rownames(tab)
     #cells <- WhichCells(object = scrna, ident = input$selectcluster)
     plot1=DimPlot(object = scrna,reduction.use=input$umapclust,group.by = input$setvar,no.legend = FALSE,do.label = TRUE, do.return=T, pt.size = input$pointclust,label.size = 7, cols.use=cpallette)
     plot2=FeaturePlot(object = scrna,reduction.use=input$umapclust, features.plot = gene, cols.use = c("grey", "blue"),do.return=T,pt.size = input$pointclust,no.legend = FALSE)
     plot2=eval(parse(text=paste("plot2$`",gene,"`",sep="")))
     plot_grid(plot1,plot2)
   })
   
   output$clustplots = renderPlot({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       clustplots()
     })
   })
   
   output$downloadclustplot <- downloadHandler(
     filename = function() {
       paste0("Clustexp_plot.pdf",sep="")
     },
     content = function(file){
       pdf(file,width=13,height = 9,useDingbats=FALSE)
       plot(clustplots())
       dev.off()
     })
   
   output$downloadclustertab <- downloadHandler(
     filename = function() { paste(input$selectcluster, '_cluster-wiseGeneExp.csv', sep='') },
     content = function(file) {
       write.csv(clustable(), file)
     })
   
   ###################################################
   ###################################################
   ################### Plot dot plot ################
   ###################################################
   ###################################################
   
   output$setdotvar = renderUI({
     scrna=fileload()
     metadata=as.data.frame(scrna@meta.data)
     metadata=metadata %>% select(starts_with("var_"))
     var=c(colnames(metadata))
     selectInput("setdotvar","Choose category",var,"pick one")
   })
   
   dotplot= reactive({
     scrna=fileload()
         validate(
           need(input$genelistfile, "Please Upload Genelist")
         )
     file=input$genelistfile
     df=fread(file$datapath,header = FALSE) #get complete gene list as string
     genes=as.vector(df$V1)
    g1=DotPlot(object = scrna, genes.plot = genes, plot.legend = TRUE,group.by=input$setdotvar,do.return=TRUE) 
    return(g1) 
   })
   output$dotplot = renderPlot({
     withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
       dotplot()
     })
   })
}#end of server