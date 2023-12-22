library(shiny)
library(Seurat)
library(SeuratObject)
library(sp)
library(shinytitle)

options(shiny.maxRequestSize = 1000 * 1024 ^ 20)

ui = fluidPage(
  use_shiny_title(),
  style = "background:#fff",
  titlePanel(title = div(img(src = "logo-sx.png", height = 100, width = 100),
                         span("Seurat Xplorer", style = "color:#527a7a"))),
  tags$style(".well {background-color:#b3b3b3;}"),
  tags$head(
    tags$style(
      type = "text/css",
      "
             #loadmessage {
               position: fixed;
               top: 0px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #A74BA9;
               z-index: 105;
             }
          "
    )
  ),
  
  conditionalPanel(condition = "$('html').hasClass('shiny-busy')",
                   tags$div("Loading...", id = "loadmessage")),
  sidebarLayout(
    sidebarPanel(
      fileInput(
        'Seurat_upload',
        'Upload Seurat file (.RDS file, max. 10 Gb)',
        accept = c('.rds', ".RDS")
      ),
      fluidRow(column(width = 2, offset = 3,
                      actionButton("uploaddata", "Upload Seurat", icon = icon("upload"),
                                   style = "background:#ccff66;font-weight:bold"))),
      hr(style = "border-color:#000; border-width: 1px;"),

      selectizeInput(inputId = "selectfeature", 
                     label = "Select Feature/Gene (not for Dim Plot)", 
                     selected = "Upload Data !", 
                     choices = "Upload Data !",
                     multiple = TRUE),
      #div(tags$em("NB: 'blend' parameter for Feature plot is set to False.")),
      
      selectizeInput("groupby", "Parameter 'group.by' (not for Feature Plot)", "Upload Data !"),
      
      selectizeInput("splitby", "Parameter 'split.by' (not for DoHeatmap)", "Upload Data !"),
      
      hr(style = "border-color:#000; border-width: 1px;"),
      fluidRow(column(width = 2, offset = 3,
                      actionButton("createplots", "Create Plots", icon = icon("chart-column"),
                                   style = "background:#99ccff;font-weight:bold"))),
      
      #hr(style = "border-color:#000; border-width: 1px;"),
      br(),
      div(style="text-align:center; color: #999999", tags$b("Copyright"),icon("copyright"),
          tags$b("2022-present"), tags$b(" | "), 
          tags$b("Altay Yuzeir"),
          tags$a(href ="https://github.com/AltayYuzeir/Shiny_Seurat_Xplorer",
                 tags$b(tags$span(style = "color: #999999", icon("github"), "GitHub")),
                 target = "_blank"), br(),
          tags$b("All rights reserved"), tags$b(" | "), tags$b("Developed with R/Shiny")),
      
      width = 4
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel(span(style = "color:#cc0000;font-weight:bold;",
                                "Violin Plot"),
                           style = "color:white",
                           plotOutput("ViolinPlot"),
                           icon = span(style = "color:#cc0000;font-weight:bold;",
                                       icon("chart-area"))
                           ),
                  tabPanel(span(style = "color:#0066ff;font-weight:bold;",
                                "Dim Plot"),
                           style = "color:black",
                           plotOutput("DimPlot"),
                           icon = span(style = "color:#0066ff;font-weight:bold;",
                                       icon("braille"))
                           ),
                  tabPanel(span(style = "color:#ff9900;font-weight:bold;",
                                "DoHeatmap"),
                           plotOutput("DoHeatMap"),
                           icon = span(style = "color:#ff9900;font-weight:bold;",
                                       icon("fire-flame-curved"))
                           ),
                  tabPanel(span(style = "color:#009933;font-weight:bold;",
                                "Feature Plot"),
                           plotOutput("FeaturePlot"),
                           icon = span(style = "color:#009933;font-weight:bold;",
                                       icon("chart-gantt"))
                           ),
      ),
     
      
    )
    
  ))

server = function(input, output, session){
  
  change_window_title(session, title = "Seurat Xplorer")
  
  seurat_upload <- reactive({
    inFile <- input$Seurat_upload
    if (is.null(inFile))
      return(NULL)
    seurat_data = readRDS(file = inFile$datapath)
    return(seurat_data)
  })
  
  
  observeEvent(input$uploaddata,{
    seurat_data = seurat_upload()
    if(!isTruthy(seurat_data)) 
      showNotification("Please provide Seurat file", type = "error")
    else {
    meta = seurat_data@meta.data
    group_by = c("--NONE--", colnames(meta))
    split_by = c("--NONE--", colnames(meta))
    
    genes = rownames(seurat_data)
    
    updateSelectizeInput(inputId = "selectfeature", choices = genes, selected = genes[1], )
    updateSelectizeInput(inputId = "groupby", choices = group_by, selected = split_by[1])
    updateSelectizeInput(inputId = "splitby", choices = split_by, selected = split_by[1])
    }
  })
  
  observeEvent(input$createplots,{
    seurat_data = seurat_upload()
    group_by = input$groupby
    split_by = input$splitby
    genes_features = input$selectfeature

    if(group_by == "--NONE--" & split_by == "--NONE--"){
      
      output$ViolinPlot = renderPlot({
        VlnPlot(object = seurat_data, features = genes_features,
                group.by = NULL, split.by = NULL, ncol = 1)
      }, height = 400*length(genes_features))
      
      output$DimPlot = renderPlot({
        DimPlot(object = seurat_data, group.by = NULL, split.by = NULL)
      })
      
      output$DoHeatMap = renderPlot({
        DoHeatmap(object = seurat_data, 
                  features = genes_features, 
                  group.by = NULL)
      }, height = 20*length(genes_features))
      
      output$FeaturePlot = renderPlot({
        FeaturePlot(object = seurat_data, features = genes_features, blend = FALSE,
                    split.by = NULL, ncol = 1)
      }, height = 400*length(genes_features))
    } else if(group_by == "--NONE--" & split_by != "--NONE--"){
      
      output$ViolinPlot = renderPlot({
        VlnPlot(object = seurat_data, features = genes_features,
                group.by = NULL, split.by = split_by, ncol = 1)
      }, height = 400*length(genes_features))
      
      output$DimPlot = renderPlot({
        DimPlot(object = seurat_data, group.by = NULL, split.by = split_by)
      })
      
      output$DoHeatMap = renderPlot({
        DoHeatmap(object = seurat_data, 
                  features = genes_features, 
                  group.by = NULL)
      }, height = 20*length(genes_features))
      
      output$FeaturePlot = renderPlot({
        FeaturePlot(object = seurat_data, features = genes_features, blend = FALSE,
                    split.by = split_by, ncol = 1)
      }, height = 400*length(genes_features))
      
    }else if(group_by != "--NONE--" & split_by == "--NONE--"){
      
      
      output$ViolinPlot = renderPlot({
        VlnPlot(object = seurat_data, features = genes_features,
                group.by = group_by, split.by = NULL, ncol = 1)
      }, height = 400*length(genes_features))
      
      output$DimPlot = renderPlot({
        DimPlot(object = seurat_data, group.by = group_by, split.by = NULL)
      })
      
      output$DoHeatMap = renderPlot({
        DoHeatmap(object = seurat_data, 
                  features = genes_features, 
                  group.by = group_by)
      }, height = 20*length(genes_features))
      
      output$FeaturePlot = renderPlot({
        FeaturePlot(object = seurat_data, features = genes_features, blend = FALSE,
                    split.by = NULL, ncol = 1)
      }, height = 400*length(genes_features))
      
    } else {
      
      output$ViolinPlot = renderPlot({
        VlnPlot(object = seurat_data, features = genes_features,
                group.by = group_by, split.by = split_by, ncol = 1)
      }, height = 400*length(genes_features))
      
      output$DimPlot = renderPlot({
        DimPlot(object = seurat_data, group.by = group_by, split.by = split_by)
      })
      
      output$DoHeatMap = renderPlot({
        DoHeatmap(object = seurat_data, 
                  features = genes_features, 
                  group.by = group_by)
      }, height = 20*length(genes_features))
      
      output$FeaturePlot = renderPlot({
        FeaturePlot(object = seurat_data, features = genes_features, blend = FALSE,
                    split.by = split_by, ncol = 1)
      }, height = 400*length(genes_features))
      
    }
  
  })
}

shinyApp(ui, server)
