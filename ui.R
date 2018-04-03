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
                     menuItem('Tsne Plot', tabName = 'tsneplot', icon = icon('hand-o-right'))
                     #          menuSubItem("PCA Plot", tabName = "dashboard"),
                     #          menuSubItem('Display Variances', tabName = 'var'),
                     #          menuSubItem('Show 3D plot', tabName = '3dplot')),
                     # menuItem('Project Summary and Results', tabName = 'summres', icon = icon('hand-o-right'), 
                     #          menuSubItem('Results', tabName = 'geneselection'),
                     #          menuSubItem('View volcano plot', tabName = 'volcanoplot'),
                     #          menuSubItem('View Limma results of Multiple Contrasts', tabName = 'multilimma'),
                     #          menuSubItem(icon=NULL,checkboxInput("check", label = "Display Contrast List", value = FALSE)))
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
    tabItem(tabName = "tsneplot",
            fluidRow(
              box(width = 8, status = "primary",title = "Tsne Plot",solidHeader = TRUE,
                  plotOutput("tsneplot",width=700,height=700)),
              box(title = "Controls",solidHeader = TRUE,width=4,status='primary',
                  selectInput("category", "Select one",c('Categories' = "var",'Cluster' = "clust", 'Gene Expression' = "geneexp"),selected = "clust"),
                  conditionalPanel(
                    condition = "input.category == 'var'",
                    uiOutput("variables")
                  ),
                  conditionalPanel(
                    condition = "input.category == 'geneexp'",
                              textInput("gene", label = h3("Gene Name"),
                                        value = "SFTPB")),#close conditional panel
                   uiOutput('range'),
                   sliderInput("pointsize", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                   #selectInput("colour", label = h3("Colour"),choices = colors(), selected = 'midnightblue'),
                   downloadButton('downloadPlot', 'Download')
                 
                    )#close control box
                    )#close fluid row
                  )#end tab item
      

      # # Second tab content
      # tabItem(tabName = "biplot",
      #         fluidRow(
      #           box(plotOutput("bigeneplot", height = 600),width=8, status='primary'),
      #           
      #           box(
      #             title = "Controls",solidHeader = TRUE,width=4,status='primary',
      #             textInput("bigene_genea", label = h3("Gene A"), 
      #                       value = "SFTPB"),
      #             #sliderInput("bigene_rangea", "Expression Limits Gene A(log2(UMI)):",
      #             #           min = 0, max = 10, value = c(0,.1),step=.01),
      #             uiOutput("bigene_rangea"),
      #             
      #             textInput("bigene_geneb", label = h3("Gene B"), 
      #                       value = "DCN"),
      #             uiOutput("bigene_rangeb"),
      #             # sliderInput("bigene_rangeb", "Expression Limits Gene B (log2(UMI)):",
      #             #            min = 0, max = 10, value = c(0,.1),step=.01),
      #             sliderInput("bigene_pointsize", "Point Size:",
      #                         min = 0, max = 5, value = 1,step=.25),
      #             radioButtons("bigene_imagetype", label = h3("Image Type"),
      #                          choices = list("PNG" = 'png', "PDF" = 'pdf'), 
      #                          selected = 'png'),
      #             downloadButton('downloadbigene', 'Download')
      #           )
      #         ),
      #         fluidRow(
      #           box(plotOutput("bigene_piechart", height = 300),width=6,status='primary'),
      #           box(tableOutput("bigene_countsTable"),
      #               width = 6, status = "primary",
      #               title = "BiGene Summary"
      #           )
      #         )#End FluidRow
      # ),#endbigeneplotTab
      # tabItem(tabName = "grpplot",
      #         fluidRow(
      #           box(plotOutput("grpPlot", height = 600),width=8),
      #           
      #           box(
      #             title = "Controls",solidHeader = TRUE,width=4,status='primary',
      #             fileInput('file1', 'Choose CSV File',
      #                       accept=c('text/csv', 
      #                                'text/comma-separated-values,text/plain', 
      #                                '.csv')),
      #             radioButtons("grpplot_imagetype", label = h3("Image Type"),
      #                          choices = list("PNG" = 'png', "PDF" = 'pdf'), 
      #                          selected = 'png'),
      #             downloadButton('downloadgrpPlot', 'Download')
      #           )
      #         )
      # ),
      # 
      # 
      # tabItem(tabName = "help",
      #         h2("Help page coming soon!")
      # )
      
    )#end of tabitems
  )#end of dashboard body
)#end of dashboard page

