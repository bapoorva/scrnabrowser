library(shinydashboard)
library(shiny)
library(shinyBS)
library(plotly)
library(d3heatmap)
library(shinyjs)
library(rglwidget)
library(SPIA)
library(cellrangerRkit)
library(reshape2)

ui <- dashboardPage(
  dashboardHeader(title = "sEuRaT",titleWidth = 350),
  dashboardSidebar(width = 350,
                   div(style="overflow-y: scroll"),
                   tags$head(tags$style(HTML(".sidebar { height: 170vh; overflow-y: auto; }" ))),
                   sidebarMenu(
                     menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
                     uiOutput("projects"),
                     menuItem('Compare tSNE Plots', tabName = 'tsneplot', icon = icon('hand-o-right')),
                     menuItem('Biplot', tabName = 'biplot', icon = icon('hand-o-right')),
                     menuItem('Differential Expression', tabName = 'deg', icon = icon('hand-o-right')),
                              #menuSubItem("Find Marker Genes", tabName = "deg")),
                     menuItem('Heatmap', tabName = 'heatmap', icon = icon('hand-o-right')),
                     menuItem('Ligand Receptor Pairs', tabName = 'ligrec2', icon = icon('hand-o-right'),
                              menuSubItem("Compare clusters", tabName = "ligrec"),
                              menuSubItem("Compare to other datasets", tabName = "compligrec"),
                              menuSubItem("Pathway Analysis", tabName = "pathway"))
                   )#end of sidebar menu
  ),#end dashboardSidebar
  
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    useShinyjs(),
    tabItems(
      tabItem(tabName = "dashboard",
              box(
                width = 10, status = "primary",solidHeader = TRUE,
                title = "scRNA Data Sets",
                tableOutput("datasetTable")
              )
      ),
    ######################################################################################################################################
    tabItem(tabName = "tsneplot",
            box(title = "Compare tSNE plots",solidHeader = TRUE,width=12,status='primary',
                fluidRow(
                  column(6,selectInput("categorya2", "Select one",c('Categories' = "var",'Cluster' = "clust", 'Gene Expression' = "geneexp"),selected = "clust")),
                  column(6,selectInput("categoryb2", "Select one",c('Categories' = "var",'Cluster' = "clust", 'Gene Expression' = "geneexp"),selected = "clust"))
                ),
                sliderInput("pointa2", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                fluidRow(
                  column(6,conditionalPanel(
                    condition = "input.categorya2 == 'var'",
                    uiOutput("tsnea2")
                  ),
                  conditionalPanel(
                    condition = "input.categorya2 == 'geneexp'",textInput("gene1a", label = "Gene Name",value = "Axin2")
                  )
                  ),
                  column(6,conditionalPanel(
                    condition = "input.categoryb2 == 'var'",
                    uiOutput("tsneb2")
                  ),
                  conditionalPanel(
                    condition = "input.categoryb2 == 'geneexp'",textInput("gene2a", label = "Gene Name",value = "Axin2")
                  )
                  )),
                fluidRow(
                  column(6,checkboxInput("subsa", label = "Check to subselect cells", value = FALSE)),
                  column(6,checkboxInput("subsb", label = "Check to subselect cells", value = FALSE))
                ),
                fluidRow(
                  column(6,conditionalPanel(
                    condition = "input.subsa ==true",uiOutput("subsaui")
                  )),
                  column(6,conditionalPanel(
                    condition = "input.subsb ==true",uiOutput("subsbui")
                  )
                  )),
                plotOutput("comptsne2", height = 600),
                downloadButton('downloadtsneplot', 'Download tSNE plot')
            )
    ),#end of degtab
    ###################################################################################################################################### 
      # tabItem(tabName = "biplot",
      #         fluidRow(
      #           box(plotOutput("bigeneplot", height = 600),width=8, status='primary',title = "Bigene Plot",solidHeader = TRUE),
      #           box(
      #             title = "Controls",solidHeader = TRUE,width=4,status='primary',
      #             textInput("bigene_genea", label = "Gene A",value = "Sox2"),
      #             textInput("bigene_geneb", label = "Gene B",value = "Sox9"),
      #             sliderInput("bigene_pointsize", "Point Size:",min = 0, max = 5, value = 1,step=.25),
      #             downloadButton('downloadbigene', 'Download bigene plot')
      #           )
      #         )#End FluidRow
      # ),#endbigeneplotTab
    
    tabItem(tabName = "biplot",
            fluidRow(
              box(plotOutput("bigeneplot", height = 600),width=8, status='primary',title = "Bigene Plot",solidHeader = TRUE),
              
              box(
                title = "Controls",solidHeader = TRUE,width=4,status='primary',
                textInput("bigene_genea", label = "Gene A",value = "Sox2"),
                #sliderInput("bigene_rangea", "Expression Limits Gene A(log2(UMI))",min = 0, max = 10, value = 0.5,step=.25),
                uiOutput("bigene_rangea"),
                textInput("bigene_geneb", label = "Gene B",value = "Sox9"),
                #sliderInput("bigene_rangeb", "Expression Limits Gene B(log2(UMI))",min = 0, max = 10, value = 1.5,step=.25),
                uiOutput("bigene_rangeb"),
                sliderInput("bigene_pointsize", "Point Size:",min = 0, max = 5, value = 1,step=.25)
                # radioButtons("bigene_imagetype", label = h3("Image Type"),
                #              choices = list("PNG" = 'png', "PDF" = 'pdf'), 
                #              selected = 'png'),
                # downloadButton('downloadbigene', 'Download')
              )
            )
    ),#endbigeneplotTab
    ######################################################################################################################################
    tabItem(tabName = "deg",
            box(title = "Compare tSNE plots",solidHeader = TRUE,width=9,status='primary',
                plotOutput("comptsne", height = 1000)
            ),

            fluidRow(
              box(
                title = "Controls",solidHeader = TRUE,width=3,status='primary',
                selectInput("categorya", "Select Category to display",c('Categories' = "var",'Cluster' = "clust"),selected = "clust"),
                conditionalPanel(
                    condition = "input.categorya == 'var'",
                    uiOutput("tsnea")
                  ),
                checkboxInput("setident", label = "Check to choose a different category to compare", value = FALSE),
                conditionalPanel(
                  condition = "input.setident ==true",uiOutput("setidentlist")
                ),
                
                uiOutput("identa"),
                uiOutput("identb"),
                sliderInput("lfc", "Log FC threshold:",min = 0.25, max = 6, value = 0.25,step=.25),
                selectInput("test", "Select test to use",c('Wilcox' = "wilcox",'T-test' = "t", 'Poisson' = "poisson",'Negative Binomial'="negbinom"),selected = "wilcox"),
                sliderInput("minpct", "Minimum Percent of cells:",min = 0.1, max = 10, value = 0.25),
                sliderInput("pointa", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                uiOutput("grptype"),
                downloadButton('downloaddeg', 'Download table'),
                downloadButton('downloadplot', 'Download Plot')
              ),
              box(DT::dataTableOutput('markergenes'),width=12, status='primary',solidHeader = TRUE,title="Marker genes")
            )#End FluidRow
    ),#end of degtab
        ######################################################################################################################################
    tabItem(tabName = "heatmap",
            box(
              title = "Controls",solidHeader = TRUE,width=12,status='primary',
              uiOutput("heatmapgenes"),
              uiOutput("hmpgrp"),
              selectInput("hmpcol", "Select one",c('PurpleYellow' = "PuYl",'BlueGreen' = "BuGn", 'RedYellow' = "RdYl", 'RedBlue'="RdBu"),selected = "geneexp"),
              downloadButton('downloadheatmap', 'Download Heatmap')
             ),
            box(plotOutput("heatmap", height = 900),width=12, status='primary',solidHeader = TRUE,title="Single cell heatmap of gene expression")
    ),#end of tab
    ######################################################################################################################################
    tabItem(tabName = "ligrec",
            box(width = 10, status = "primary",solidHeader = TRUE,title = "Controls",
                radioButtons("clust","Select one", c("All clusters"="all","Select Cluster"="clust"),selected = "clust"),
                radioButtons("gene","Select one", c("All genes"="allgene","Enter Genelist"="genelist"),selected = "allgene"),
                
                conditionalPanel(
                  condition = "input.clust == 'all' && input.gene == 'genelist'" ,
                  uiOutput("list1"),
                  uiOutput("list2")
                ),
                conditionalPanel(
                  condition = "input.clust == 'clust' && input.gene == 'allgene'" ,
                  uiOutput("clust1"),
                  uiOutput("clust2")
                ),
                conditionalPanel(
                  condition = "input.clust == 'clust' && input.gene == 'genelist'" ,
                  uiOutput("list1.1"),
                  uiOutput("list2.1"),
                  uiOutput("clust1.1"),
                  uiOutput("clust2.1")
                )
            ),
            
            box(
              width = 10, status = "primary",solidHeader = TRUE,
              title = "Ligand Receptor pairs",
              DT::dataTableOutput('pairs_res')
            )#end of box
    ),#end of tabitem
    ######################################################################
    ######################################################################
    tabItem(tabName = "compligrec",
            box(width = 6, status = "primary",solidHeader = TRUE,title = "Ligand Selection Panel",
                selectInput("ligand", "Select Experiment type",c('RNA-Seq' = "rna",'Single Cell' = "scrna", 'Microarray' = "microarray")),
                uiOutput("ligprj"),
                uiOutput("ligtype"),
                conditionalPanel(
                  condition = "input.ligand == 'rna' | input.ligand == 'microarray'" ,
                  sliderInput("explig", label = "Set Expression threshold", min =6,max = 12, value = 6)),
                conditionalPanel(
                  condition = "input.ligand == 'scrna'" ,
                  sliderInput("ligumi", label = "Set UMI threshold", min =1,max = 25, value = 1),
                  sliderInput("ligsamp", label = "Set Percent Samples", min =0,max = 100, value = 50)),
                checkboxInput("liggene", label = "Upload Gene List", value = FALSE),
                
                conditionalPanel(
                  condition = "input.liggene ==true",
                  fileInput('liggeneli', 'Upload Receptor Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
                )),
            box(width = 6, status = "primary",solidHeader = TRUE,title = "Receptor Selection Panel",
                selectInput("receptor", "Select Experiment type",c('RNA-Seq' = "rna",'Single Cell' = "scrna", 'Microarray' = "microarray")),
                uiOutput("recprj"),
                uiOutput("rectype"),
                conditionalPanel(
                  condition = "input.receptor == 'rna' | input.receptor == 'microarray'" ,
                  sliderInput("exprec", label = "Set Expression threshold", min =6,max = 12, value = 6)),
                conditionalPanel(
                  condition = "input.receptor == 'scrna'" ,
                  sliderInput("recumi", label = "Set UMI threshold", min =1,max = 25, value = 1),
                  sliderInput("recsamp", label = "Set Percent Samples", min =0,max = 100, value = 50)),
                checkboxInput("recgene", label = "Upload Gene List", value = FALSE),
                conditionalPanel(
                  condition = "input.recgene ==true",
                  fileInput('recgeneli', 'Upload Receptor Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
                )),
            box(width = 12, status = "primary",solidHeader = TRUE,title = "Ligand-Receptor pairs",
                DT::dataTableOutput('ligrecpairs'),uiOutput("dwldtab"))
    ),#end of tabitem
    ######################################################################
    ######################################################################
    tabItem(tabName = "pathway",
            box(width = 12, status = "primary",solidHeader = TRUE,title = "Ligand Receptor Pairs",
                DT::dataTableOutput('rec')
            ),#end of box
            box(width = 12, status = "primary",solidHeader = TRUE,title = "KEGG Pathways",
                DT::dataTableOutput('Keggpaths')
            ),
            box(width = 12,height = 12, status = "primary",solidHeader = TRUE,title = "Pathview",
                plotOutput("plots")
            )
    )#end of tabitem
    ######################################################################
    ######################################################################
    )#end of tabitems
  )#end of dashboard body
)#end of dashboard page

