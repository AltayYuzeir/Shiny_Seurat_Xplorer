library(shiny)
library(Seurat)
library(SeuratObject)
library(sp)

options(shiny.maxRequestSize = 1000 * 1024 ^ 20)

ui = fluidPage(
  titlePanel(title = div(img(src = "logo-sx.png", height = 100, width = 100),
                         "Seurat Xplorer")),
  sidebarLayout(
    sidebarPanel(
      fileInput(
        'Seurat_upload',
        'Upload Seurat file (.RDS file, max. 10 Gb)',
        accept = c('.rds', ".RDS")
      ),
      fluidRow(column(width = 2, offset = 3,
                      actionButton("uploaddata", "Upload Seurat", icon = icon("upload"),
                                   style = "background:#ccff66;"))),
      hr(),
      div(tags$em("NB: Max of 9 features can be selected")),
      
      selectizeInput(inputId = "selectfeature", 
                     label = "Select Feature/Gene (not for Dim Plot)", 
                     selected = "Upload Data !", 
                     choices = "Upload Data !",
                     multiple = TRUE),
      
      selectizeInput("groupby", "Parameter 'group.by' (not for Feature Plot)", "Upload Data !"),
      div(tags$em("NB: 'blend' parameter for Feature plot is set to False.")),
      
      selectizeInput("splitby", "Parameter 'split.by' (not for DoHeatmap)", "Upload Data !"),
      
      hr(),
      fluidRow(column(width = 2, offset = 3,
                      actionButton("createplots", "Create Plots", icon = icon("chart-column"),
                                   style = "background:#99ccff;"))),
      
      tags$head(tags$style(
        HTML("hr {border-top: 1px solid #000000;}")
      )),
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
      
      width = 4
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Violin Plot", plotOutput("ViolinPlot")),
                  tabPanel("Dim Plot", plotOutput("DimPlot")),
                  tabPanel("DoHeatmap", plotOutput("DoHeatMap")),
                  tabPanel("Feature Plot", plotOutput("FeaturePlot")),
      ),
      hr(),
      div(style = "text-align:center; color: #6699ff", 
          "Developed by Altay Yuzeir, Bonn-Germany, 2022",
          tags$a(href ="https://github.com/AltayYuzeir/Seurat_Xplorer/",
                 tags$b(tags$span(style = "color: #6699ff", icon("github"), "GitHub")),
                 target = "_blank"))
      
    )
    
  ))

server = function(input, output, session){
  
  seurat_upload <- reactive({
    inFile <- input$Seurat_upload
    if (is.null(inFile))
      return(NULL)
    seurat_data = readRDS(file = inFile$datapath)
    return(seurat_data)
  })
  
  
  observeEvent(input$uploaddata,{
    seurat_data = seurat_upload()
    meta = seurat_data@meta.data
    group_by = c("--NONE--", colnames(meta))
    split_by = c("--NONE--", colnames(meta))
    
    genes = rownames(seurat_data)
    
    updateSelectizeInput(inputId = "selectfeature", choices = genes, selected = "IGKC", )
    updateSelectizeInput(inputId = "groupby", choices = group_by, selected = "Time")
    updateSelectizeInput(inputId = "splitby", choices = split_by, selected = split_by[1])
    
  })
  
  observeEvent(input$createplots,{
    seurat_data = seurat_upload()
    group_by = input$groupby
    split_by = input$splitby
    genes_features = input$selectfeature
    if(length(genes_features) > 9) genes_features = genes_features[1:9]
    
    if(group_by == "--NONE--" & split_by == "--NONE--"){
      
      output$ViolinPlot = renderPlot({
        VlnPlot(object = seurat_data, features = genes_features,
                group.by = NULL, split.by = NULL)
      })
      
      output$DimPlot = renderPlot({
        DimPlot(object = seurat_data, group.by = NULL, split.by = NULL)
      })
      
      output$DoHeatMap = renderPlot({
        DoHeatmap(object = seurat_data, 
                  features = genes_features, 
                  group.by = NULL)
      })
      
      output$FeaturePlot = renderPlot({
        FeaturePlot(object = seurat_data, features = genes_features, blend = FALSE,
                    split.by = NULL)
      })
    } else if(group_by == "--NONE--" & split_by != "--NONE--"){
      
      output$ViolinPlot = renderPlot({
        VlnPlot(object = seurat_data, features = genes_features,
                group.by = NULL, split.by = split_by)
      })
      
      output$DimPlot = renderPlot({
        DimPlot(object = seurat_data, group.by = NULL, split.by = split_by)
      })
      
      output$DoHeatMap = renderPlot({
        DoHeatmap(object = seurat_data, 
                  features = genes_features, 
                  group.by = NULL)
      })
      
      output$FeaturePlot = renderPlot({
        FeaturePlot(object = seurat_data, features = genes_features, blend = FALSE,
                    split.by = split_by)
      })
      
    }else if(group_by != "--NONE--" & split_by == "--NONE--"){
      
      
      output$ViolinPlot = renderPlot({
        VlnPlot(object = seurat_data, features = genes_features,
                group.by = group_by, split.by = NULL)
      })
      
      output$DimPlot = renderPlot({
        DimPlot(object = seurat_data, group.by = group_by, split.by = NULL)
      })
      
      output$DoHeatMap = renderPlot({
        DoHeatmap(object = seurat_data, 
                  features = genes_features, 
                  group.by = group_by)
      })
      
      output$FeaturePlot = renderPlot({
        FeaturePlot(object = seurat_data, features = genes_features, blend = FALSE,
                    split.by = NULL)
      })
      
    } else {
      
      output$ViolinPlot = renderPlot({
        VlnPlot(object = seurat_data, features = genes_features,
                group.by = group_by, split.by = split_by)
      })
      
      output$DimPlot = renderPlot({
        DimPlot(object = seurat_data, group.by = group_by, split.by = split_by)
      })
      
      output$DoHeatMap = renderPlot({
        DoHeatmap(object = seurat_data, 
                  features = genes_features, 
                  group.by = group_by)
      })
      
      output$FeaturePlot = renderPlot({
        FeaturePlot(object = seurat_data, features = genes_features, blend = FALSE,
                    split.by = split_by)
      })
      
    }
    
  })
}

shinyApp(ui, server)
