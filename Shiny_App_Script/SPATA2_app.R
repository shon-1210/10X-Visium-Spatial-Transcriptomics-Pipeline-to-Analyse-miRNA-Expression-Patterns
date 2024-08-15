library(shiny)
library(shinyjs)
library(BiocManager)
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(dplyr)
library(SPATAData)
library(SPATA2)
library(stringr)
library(readr)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)

# Define the directory containing .rds files
rds_directory <- "data/RDS_Files_V3_Script/"

# Define the path to the target file
targets_file_path <- "data/Targetscan_Data/MMS_Predicted_Targets_Context_Scores.default_predictions.txt"

# Function to list .rds files in the directory and its subdirectories
list_rds_files <- function(directory) {
  files <- list.files(directory, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
  file_names <- basename(files)
  names(files) <- file_names
  return(files)
}

ui <- fluidPage(
  useShinyjs(),  # Initialize shinyjs
  tags$head(
    tags$style(HTML("
      body {
        zoom: 75%;
        -webkit-transform: scale(1);
        transform: scale(1);
        transform-origin: 0 0;
      }
      
      @media only screen and (max-width: 768px) {
        body {
          zoom: 100%;
          -webkit-transform: scale(1);
          transform: scale(1);
        }
      }
      hr.custom-line {
        border: 1px solid black; /* Set the horizontal line color to black */
      }
      .custom-text {
        color: #505050; /* Darker gray color for descriptive text */
      }
    ")),
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1.0, user-scalable=no")
  ),
  titlePanel("SPATA2 Dataset Specific Analysis Heatmap Generator"),
  sidebarLayout(
    sidebarPanel(
      helpText("Please wait while the webpage loads. Then, select the desired .RDS file to generate the heatmaps."),
      selectInput("spataFile", "Choose SPATA2 RDS File", choices = list_rds_files(rds_directory), selected = ""),
      selectInput("clusteringMethod", "Clustering Method", choices = c("K-means", "BayesSpace")),
      sliderInput("topN", "Number of Top Targets (topN)", min = 100, max = 1000, value = 100, step = 100),
      actionButton("generateHeatmap", "Generate Heatmap"),
      tags$hr(class = "custom-line"),  # Dividing line with custom color
      helpText("To generate the comparison heatmap, please select pairs of clusters for comparison from the dropdown menus below after the first heatmap has been generated. Scroll down the webpage to view this heatmap."),
      fluidRow(
        column(6, selectInput("cluster1_1", "Select Cluster 1 (Pair 1)", choices = NULL)),
        column(6, selectInput("cluster2_1", "Select Cluster 2 (Pair 1)", choices = NULL))
      ),
      fluidRow(
        column(6, selectInput("cluster1_2", "Select Cluster 1 (Pair 2)", choices = NULL)),
        column(6, selectInput("cluster2_2", "Select Cluster 2 (Pair 2)", choices = NULL))
      ),
      fluidRow(
        column(6, selectInput("cluster1_3", "Select Cluster 1 (Pair 3)", choices = NULL)),
        column(6, selectInput("cluster2_3", "Select Cluster 2 (Pair 3)", choices = NULL))
      ),
      actionButton("generateComparisonHeatmap", "Generate Comparison Heatmap"),
      helpText("Warning: Changing the zoom level of the current browser tab will result in the heatmaps being resized to fit the chosen zoom level.")
    ),
    mainPanel(
      plotOutput("combinedPlot"),
      plotOutput("heatmapPlot"),
      plotOutput("comparisonHeatmap"),
      verbatimTextOutput("debug")
    )
  )
)

server <- function(input, output, session) {
  # Disable inputs initially
  disable("clusteringMethod")
  disable("topN")
  disable("generateHeatmap")
  disable("cluster1_1")
  disable("cluster2_1")
  disable("cluster1_2")
  disable("cluster2_2")
  disable("cluster1_3")
  disable("cluster2_3")
  disable("generateComparisonHeatmap")
  
  observe({
    rds_files <- list_rds_files(rds_directory)
    updateSelectInput(session, "spataFile", choices = rds_files)
  })
  
  observeEvent(input$spataFile, {
    # Enable inputs when an RDS file is selected
    if (!is.null(input$spataFile) && input$spataFile != "") {
      enable("clusteringMethod")
      enable("topN")
      enable("generateHeatmap")
    }
  })
  
  data <- reactive({
    req(input$spataFile, input$clusteringMethod)
    
    spataObj <- readRDS(input$spataFile)
    sampleName <- getSampleName(spataObj)
    
    # Get and log1p normalize the count matrix
    exp_matrix_scaled <- getCountMatrix(spataObj)
    exp_matrix_scaled <- as.data.frame(log1p(exp_matrix_scaled))
    
    featureDf <- getFeatureDf(spataObj)
    
    if (input$clusteringMethod == "K-means") {
      cluTab <- featureDf %>%
        dplyr::select(barcodes, kmeans_4_HW) %>%
        column_to_rownames(var = "barcodes")
      cluTab$clusters <- paste0("cluster.", cluTab$kmeans_4_HW)
    } else {
      cluTab <- featureDf %>%
        dplyr::select(barcodes, bayes_space) %>%
        column_to_rownames(var = "barcodes")
      cluTab$clusters <- paste0("cluster.", cluTab$bayes_space)
    }
    
    tarTab <- read.table(targets_file_path, sep = "\t", header = TRUE, check.names = FALSE)
    tarTab <- tarTab[tarTab$`Gene Symbol` %in% rownames(exp_matrix_scaled), ]
    tarTab <- tarTab[grepl("^mmu-", tarTab$miRNA), ]
    tarTab <- tarTab[order(tarTab$`weighted context++ score`, tarTab$`context++ score`), ]
    
    # Enable cluster selection inputs after data is ready
    enable("cluster1_1")
    enable("cluster2_1")
    enable("cluster1_2")
    enable("cluster2_2")
    enable("cluster1_3")
    enable("cluster2_3")
    enable("generateComparisonHeatmap")
    
    list(expTab = exp_matrix_scaled, cluTab = cluTab, tarTab = tarTab, sampleName = sampleName)
  })
  
  observe({
    req(data())
    clusters <- sort(unique(data()$cluTab$clusters))
    updateSelectInput(session, "cluster1_1", choices = clusters)
    updateSelectInput(session, "cluster2_1", choices = clusters)
    updateSelectInput(session, "cluster1_2", choices = clusters)
    updateSelectInput(session, "cluster2_2", choices = clusters)
    updateSelectInput(session, "cluster1_3", choices = clusters)
    updateSelectInput(session, "cluster2_3", choices = clusters)
  })
  
  output$combinedPlot <- renderPlot({
    req(data())
    spataObj <- readRDS(input$spataFile)
    
    plot1 <- plotSurface(
      object = spataObj,
      color_by = "bayes_space",
      pt_size = 2.3,
      pt_clrp = "uc"
    ) + 
      labs(subtitle = "BayesSpace Clusters") +
      theme(
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)
      )
    
    plot2 <- plotSurface(
      object = spataObj,
      color_by = "kmeans_4_HW",
      pt_clrp = "uc",
      pt_size = 2.3
    ) + 
      labs(subtitle = "HW K-means Clusters") +
      theme(
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)
      )
    
    combined_plot <- plot1 | plot2
    print(combined_plot)
  })
  
  output$heatmapPlot <- renderPlot({
    req(data())
    expTab <- data()$expTab
    cluTab <- data()$cluTab
    tarTab <- data()$tarTab
    
    mirs <- c("mmu-miR-124-3p.1", "mmu-miR-9-5p", "mmu-miR-122-5p", "mmu-miR-1a-3p", "mmu-miR-133a-3p.1",
              "mmu-let-7a-5p")
    
    meanExpr <- data.frame(row.names = rownames(expTab))
    for (cluster in sort(unique(cluTab$clusters))) {
      barcodes <- rownames(cluTab)[cluTab$clusters == cluster]
      if (length(barcodes) == 0) next
      meanExpr[, cluster] <- rowMeans(expTab[, barcodes, drop = FALSE])
    }
    
    foldChange <- data.frame(row.names = rownames(meanExpr))
    for (cluster in colnames(meanExpr)) {
      foldChange[, cluster] <- (meanExpr[, cluster]) - (rowMeans(meanExpr[, colnames(meanExpr) != cluster, drop = FALSE]))
    }
    
    pvalTab <- data.frame(row.names = mirs)
    
    for (mir in mirs) {
      tarGenesAll <- unique(tarTab$`Gene Symbol`[tarTab$miRNA == mir])
      tarGenesTop <- tarGenesAll[seq_len(min(input$topN, length(tarGenesAll)))]
      for (cluster in colnames(foldChange)) {
        x <- foldChange[tarGenesTop, cluster]
        y <- foldChange[!rownames(foldChange) %in% tarGenesTop, cluster]
        pval <- wilcox.test(x, y)$p.value
        sign_TarvsNonTar <- sign(median(x) - median(y))
        pval <- log10(pval) * -sign_TarvsNonTar
        pvalTab[mir, cluster] <- signif(pval, 4)
      }
    }
    
    logpvalMatrix <- as.matrix.data.frame(pvalTab)
    logpvalMatrix[logpvalMatrix > 10] <- 10
    logpvalMatrix[logpvalMatrix < -10] <- -10
    
    breaksList <- seq(-10, 10, by = 0.1)
    myColorPalette <- colorRampPalette(c("darkolivegreen3", "white", "red"))(length(breaksList) - 1)
    
    pheatmap(logpvalMatrix, cluster_rows = FALSE, cluster_cols = FALSE,
             display_numbers = TRUE,
             main = "Heatmap of transformed log10(p-values) for miRNA targets across clusters",
             color = myColorPalette,
             breaks = breaksList,
             fontsize_number = 10,
             angle_col = 45)
  })
  
  output$comparisonHeatmap <- renderPlot({
    req(data())
    expTab <- data()$expTab
    cluTab <- data()$cluTab
    tarTab <- data()$tarTab
    
    mirs <- c("mmu-miR-124-3p.1", "mmu-miR-9-5p", "mmu-miR-122-5p", "mmu-miR-1a-3p", "mmu-miR-133a-3p.1",
              "mmu-let-7a-5p")
    
    meanExpr <- data.frame(row.names = rownames(expTab))
    for (cluster in sort(unique(cluTab$clusters))) {
      barcodes <- rownames(cluTab)[cluTab$clusters == cluster]
      if (length(barcodes) == 0) next
      meanExpr[, cluster] <- rowMeans(expTab[, barcodes, drop = FALSE])
    }
    
    foldChange <- data.frame(row.names = rownames(meanExpr))
    cluster_pairs <- list(
      c(input$cluster1_1, input$cluster2_1),
      c(input$cluster1_2, input$cluster2_2),
      c(input$cluster1_3, input$cluster2_3)
    )
    
    pvalTab <- data.frame(row.names = mirs)
    
    for (specific_clusters in cluster_pairs) {
      if (!all(specific_clusters %in% colnames(meanExpr))) {
        stop("One or more specified clusters not found in the mean expression data.")
      }
      
      foldChange[[paste(specific_clusters[1], "vs", specific_clusters[2])]] <- meanExpr[, specific_clusters[1]] - meanExpr[, specific_clusters[2]]
      
      for (mir in mirs) {
        tarGenesAll <- unique(tarTab$`Gene Symbol`[tarTab$miRNA == mir])
        tarGenesTop <- tarGenesAll[seq_len(min(input$topN, length(tarGenesAll)))]
        
        x <- foldChange[tarGenesTop, paste(specific_clusters[1], "vs", specific_clusters[2])]
        y <- foldChange[!rownames(foldChange) %in% tarGenesTop, paste(specific_clusters[1], "vs", specific_clusters[2])]
        pval <- wilcox.test(x, y)$p.value
        sign_TarvsNonTar <- sign(median(x) - median(y))
        pval <- log10(pval) * -sign_TarvsNonTar
        pvalTab[mir, paste(specific_clusters[1], "vs", specific_clusters[2])] <- signif(pval, 4)
      }
    }
    
    logpvalMatrix <- as.matrix.data.frame(pvalTab)
    
    breaksList <- seq(-10, 10, by = 0.1)
    myColorPalette <- colorRampPalette(c("darkolivegreen3", "white", "red"))(length(breaksList) - 1)
    
    pheatmap(logpvalMatrix, cluster_rows = FALSE, cluster_cols = FALSE,
             display_numbers = TRUE,
             main = "Heatmap of transformed log10(p-values) for miRNA targets between clusters",
             color = myColorPalette,
             breaks = breaksList,
             fontsize_number = 10,
             angle_col = 45)
  })
  
  output$debug <- renderPrint({
    str(data())
  })
}

shinyApp(ui = ui, server = server)
