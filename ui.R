library(shinydashboard)
#library(shinyIncubator)
library(shiny)
library(shinyBS)
library(plotly)
library(d3heatmap)
library(shinyjs)
library(rglwidget)
library(SPIA)
library(cellrangerRkit)
library(reshape2)
source('plots.R')
#load('data/jarod_scRNA.rdata')


ui <- dashboardPage(
  dashboardHeader(title = "sEuRaT",titleWidth = 350),
  dashboardSidebar(width = 350,
                   div(style="overflow-y: scroll"),
                   tags$head(tags$style(HTML(".sidebar { height: 170vh; overflow-y: auto; }" ))),
                   sidebarMenu(
                     menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
                     uiOutput("projects"),
                     menuItem('tSNE Plot', tabName = 'tsneplot', icon = icon('hand-o-right')),
                     menuItem('Biplot', tabName = 'biplot', icon = icon('hand-o-right')),
                     menuItem('Differential Expression', tabName = 'mgenes', icon = icon('hand-o-right'),
                              menuSubItem("Find Marker Genes", tabName = "deg"),
                              menuSubItem("Violin Plots", tabName = "violinplot"),
                              menuSubItem("Feature plots", tabName = "featureplots")),
                     #menuItem('Violin Plots', tabName = 'violinplot', icon = icon('hand-o-right')),
                     menuItem('Heatmap', tabName = 'heatmap', icon = icon('hand-o-right'))
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
            fluidRow(
              box(width = 8, status = "primary",title = "Tsne Plot",solidHeader = TRUE,
                  plotlyOutput("tsneplot",height=700)),
              box(title = "Controls",solidHeader = TRUE,width=4,status='primary',
                  selectInput("category", "Select one",c('Categories' = "var",'Cluster' = "clust", 'Gene Expression' = "geneexp"),selected = "geneexp"),
                  conditionalPanel(
                    condition = "input.category == 'var'",
                    uiOutput("variables")
                  ),
                  conditionalPanel(
                    condition = "input.category == 'geneexp'",textInput("gene", label = "Gene Name",value = "Axin2")
                    ),#close conditional panel
                   uiOutput('range'),
                   downloadButton('downloadPlot', 'Download')
                    )#close control box
                    )#close fluid row
                  ),#end tab item
    ###################################################################################################################################### 
      tabItem(tabName = "biplot",
              fluidRow(
                box(plotOutput("bigeneplot", height = 600),width=8, status='primary'),
                box(
                  title = "Controls",solidHeader = TRUE,width=4,status='primary',
                  textInput("bigene_genea", label = "Gene A",value = "Sox2"),
                  textInput("bigene_geneb", label = "Gene B",value = "Sox9"),
                  sliderInput("bigene_pointsize", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                  downloadButton('downloadbigene', 'Download')
                )
              )#End FluidRow
      ),#endbigeneplotTab
    ######################################################################################################################################
    tabItem(tabName = "deg",
            box(title = "Compare tSNE plots",solidHeader = TRUE,width=12,status='primary',
                fluidRow(
                  column(6,selectInput("categorya", "Select one",c('Categories' = "var",'Cluster' = "clust", 'Gene Expression' = "geneexp"),selected = "clust")),
                  column(6,selectInput("categoryb", "Select one",c('Categories' = "var",'Cluster' = "clust", 'Gene Expression' = "geneexp"),selected = "clust"))
                  # column(6,uiOutput("tsnea")),
                  # column(6,uiOutput("tsneb"))
                ),
                fluidRow(
                  column(6,sliderInput("pointa", "Point Size:",min = 0, max = 5, value = 1,step=.25)),
                  column(6,sliderInput("pointb", "Point Size:",min = 0, max = 5, value = 1,step=.25))
                ),
                fluidRow(
                  column(6,conditionalPanel(
                    condition = "input.categorya == 'var'",
                    uiOutput("tsnea")
                  ),
                  conditionalPanel(
                    condition = "input.categorya == 'geneexp'",textInput("gene1", label = "Gene Name",value = "Axin2")
                  )
                  ),
                  column(6,conditionalPanel(
                    condition = "input.categoryb == 'var'",
                    uiOutput("tsneb")
                  ),
                  conditionalPanel(
                    condition = "input.categoryb == 'geneexp'",textInput("gene2", label = "Gene Name",value = "Axin2")
                  )
                  )),
                plotOutput("comptsne", height = 900)
            ),
            fluidRow(
              box(
                title = "Controls",solidHeader = TRUE,width=4,status='primary',
                textInput("identa", label = "Identity A",value = "1"),
                textInput("identb", label = "Identity B"),
                sliderInput("lfc", "Log FC threshold:",min = 0.25, max = 6, value = 0.25,step=.25),
                selectInput("test", "Select test to use",c('Wilcox' = "wilcox",'T-test' = "t", 'Poisson' = "poisson",'Negative Binomial'="negbinom","DESeq2"="DESeq2"),selected = "wilcox"),
                sliderInput("minpct", "Minimum Percent of cells:",min = 0.1, max = 10, value = 0.25),
                downloadButton('downloaddeg', 'Download table')
              ),
              box(DT::dataTableOutput('markergenes'),width=8, status='primary',solidHeader = TRUE,title="Marker genes")
            )#End FluidRow
    ),#end of degtab
    ######################################################################################################################################
    tabItem(tabName = "violinplot",
              box(
                title = "Controls",solidHeader = TRUE,width=12,status='primary',
                sliderInput("vplot", "Number of top genes to plot:",min = 1, max = 20, value = 4),
                uiOutput("grptype"),
                downloadButton('downloadviolin', 'Download')
              ),
            box(plotOutput("violinplot", height = 900),width=12, status='primary',solidHeader = TRUE,title="Top genes- Violin Plot")
    ),#end of tab
    ######################################################################################################################################
    tabItem(tabName = "featureplots",
            box(
              title = "Controls",solidHeader = TRUE,width=12,status='primary',
              sliderInput("cowplot", "Number of top genes to plot:",min = 1, max = 16, value = 4),
              sliderInput("marker_pointsize", "Point Size:",min = 0, max = 5, value = 1,step=.25),
              downloadButton('dwldfeature', 'Download plot')
            ),
            box(plotOutput("markgeneplot", height = 900),width=12, status='primary',solidHeader = TRUE,title="Top Marker genes")
    ),#end of tab
    ######################################################################################################################################
    tabItem(tabName = "heatmap",
            box(
              title = "Controls",solidHeader = TRUE,width=12,status='primary',
              uiOutput("heatmapgenes"),
              uiOutput("hmpgrp"),
              selectInput("hmpcol", "Select one",c('PurpleYellow' = "PuYl",'BlueGreen' = "BuGn", 'RedYellow' = "RdYl", 'RedBlue'="RdBu"),selected = "geneexp")
             ),
            box(plotOutput("heatmap", height = 900),width=12, status='primary',solidHeader = TRUE,title="Single cell heatmap of gene expression")
    )#end of tab
    ######################################################################################################################################
    )#end of tabitems
  )#end of dashboard body
)#end of dashboard page

