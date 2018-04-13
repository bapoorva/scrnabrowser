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
                     menuItem('Differential Expression', tabName = 'mgenes', icon = icon('hand-o-right'),
                              menuSubItem("Find Marker Genes", tabName = "deg")),
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
      tabItem(tabName = "biplot",
              fluidRow(
                box(plotOutput("bigeneplot", height = 600),width=8, status='primary',title = "Bigene Plot",solidHeader = TRUE),
                box(
                  title = "Controls",solidHeader = TRUE,width=4,status='primary',
                  textInput("bigene_genea", label = "Gene A",value = "Sox2"),
                  textInput("bigene_geneb", label = "Gene B",value = "Sox9"),
                  sliderInput("bigene_pointsize", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                  downloadButton('downloadbigene', 'Download bigene plot')
                )
              )#End FluidRow
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
                sliderInput("pointa", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                textInput("identa", label = "First cluster to compare",value = "0"),
                textInput("identb", label = "Second Cluster(s) to compare", value = "1"),
                sliderInput("lfc", "Log FC threshold:",min = 0.25, max = 6, value = 0.25,step=.25),
                selectInput("test", "Select test to use",c('Wilcox' = "wilcox",'T-test' = "t", 'Poisson' = "poisson",'Negative Binomial'="negbinom","DESeq2"="DESeq2"),selected = "wilcox"),
                sliderInput("minpct", "Minimum Percent of cells:",min = 0.1, max = 10, value = 0.25),
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
    )#end of tab
    ######################################################################################################################################
    )#end of tabitems
  )#end of dashboard body
)#end of dashboard page

