# 10X Visium Spatial Transcriptomics Pipeline to Analyse miRNA Expression Patterns

**Abstract:**   
MicroRNAs (miRNAs) regulate messenger RNA (mRNA) expression post-transcriptionally in mammals by non-slicing repression. When a miRNA forms the RNA-induced silencing complex (RISC) with an Argonaute (AGO) protein, it targets specific mRNAs 
which are known as its targets. Since their discovery in 1993 by Lee et al. in C. elegans, the functions of miRNAs have been widely debated, with various mechanisms proposed. Stark et al. found in Drosophila that certain miRNAs and their targets are not co-expressed but are instead found in separate cells, suggesting a regulatory decision at specific transition points. This project investigates whether miRNAs act as guardians of gene expression between neighbouring cell types by analysing the expression levels of their targets at these transition points.  
To explore this, we developed a pipeline to detect the effect of tissue-specific miRNAs in Spatial Transcriptomic (ST) datasets, where gene expression levels are mapped to tissue locations. Selected ST datasets were then analysed to find neighbouring cell types with high and low expression of miRNA targets, indicating potential miRNA activity. By visualising the spatial expression patterns of the targets at these transition points using ST analysis tools, evidence can be collected to support the hypothesis. miR-124 and miR-1, specifically expressed in Brain and Heart tissue respectively, were detected with high confidence and with neighbouring clusters of high and low expression obtained for both. However, the hypothesis could not be validated in the brain dataset that was analysed in detail, as the expression of the targets was low only in one region and expressed to a higher level in varying degrees for all other cellular regions. The results however were promising as an initial step and should be followed up on by testing more datasets.

## Key Features

- **End-to-end spatial analysis pipeline:** Process raw 10x Visium data (gene-barcode matrices and tissue images) into analysis-ready form using Seurat and SPATA2 frameworks.  
- **Denoising with autoencoders:** Improve data quality by applying an autoencoder neural network to remove technical noise and outliers from the gene expression matrix.  
- **Flexible clustering methods:** Identify spatial domains using both conventional K-means clustering and spatially aware clustering with BayesSpace, to capture cell-type clusters and tissue regions.  
- **miRNA target integration:** Incorporate TargetScan predicted miRNA–mRNA interactions to examine if differential gene expression at cluster boundaries could be explained by miRNA regulation.  
- **Spatial expression visualizations:** Generate heatmaps and spatial plots overlaying gene expression or cluster information on the tissue image, highlighting patterns (e.g. cluster-specific marker expression, cluster vs. cluster comparisons).  
- **Interactive Shiny app:** Includes a Shiny web application for interactive analysis – allowing users to select clusters or genes, visualize spatial expression maps, and explore miRNA target effects in the tissue context.  

## How It Works

1. **Data Input & Setup:** The pipeline starts by loading a spatial transcriptomics dataset (10x Visium output, including the count matrix and spatial coordinates/image) and initializing it as a SPATA2/Seurat object. Basic preprocessing like quality control and normalization (using Seurat’s SCTransform v2) is applied to prepare the data for analysis.  
2. **Denoising with Autoencoders:** An autoencoder neural network is then trained on the expression data to learn a compressed representation. The network’s output is used to reconstruct a denoised expression matrix, reducing technical noise and emphasizing true biological signal.  
3. **Clustering (K-means & BayesSpace):** The pipeline performs unsupervised clustering on the spatial spots in two ways: (a) **K-means clustering** on the gene expression features to group spots by transcriptional similarity, and (b) **BayesSpace clustering**, which leverages a Bayesian model and spatial information to group neighboring spots into the same cluster. These approaches help identify distinct cell-type regions or spatial domains within the tissue.  
4. **miRNA Target Analysis:** For each identified cluster and boundary region, the pipeline analyzes gene expression differences and cross-references them with known miRNA target data. Using TargetScan predictions, it checks whether the genes up- or down-regulated at cluster edges are enriched for targets of certain miRNAs (especially tissue-specific miRNAs). This helps evaluate the hypothesis that those miRNAs might be regulating genes at the interface of different cell types.  
5. **Visualization & Outputs:** The results are compiled into easily interpretable outputs. The pipeline generates tables of cluster-specific marker genes and miRNA target overlaps, as well as **spatial heatmaps** that illustrate where particular gene expression or miRNA activity is high or low across the tissue slice. For example, it can produce a heatmap showing a cluster’s top genes compared to all other clusters, or pairwise cluster comparison plots highlighting boundary effects. These visualizations are output as HTML reports and/or image files for further inspection.  
6. **Interactive Exploration:** Finally, the included Shiny app allows users to interactively explore the processed data. Users can load a processed SPATA2 object (as an `.rds` file), then dynamically generate plots and compare clusters via a web interface – without needing to rerun code. This is especially useful for examining specific genes or miRNAs of interest in the spatial context, beyond the static results.  

## Pipeline Flowchart

Below is a high-level flowchart of the pipeline:

![Pipeline Flowchart](figures/pipeline_flowchart.svg)


```mermaid
flowchart TB
  %% Compact vertical pipeline (GitHub-safe)

  %% Nodes
  start((Start))
  input[Input\n10x Visium outputs\nfiltered_feature_bc_matrix.h5 + spatial/*]
  validate[Validate inputs\n- paths & image\n- gene-barcode matrix]
  reqOK{Required files present?}
  fixInputs[Fix inputs\nthen re-validate]

  init[Init SPATA2 / Seurat\nSCTransform v2]
  ae[Autoencoder assess\nchoose activation/bottleneck]
  denoise[Denoise matrix\nmake denoised layer]

  cluster[Clustering (both)\nBayesSpace + K-means (HW)]
  freeze[Freeze results\nsave .RDS]
  enoughK{>= 3 clusters?}
  tuneClust[Tune clustering\nadjust BayesSpace q / params\nre-run K-means]

  targets[Targets (TargetScan)\nchoose topN; add let-7 control]
  logfc[Compute logFC\ncluster vs rest; pairwise; distant ctrl]
  wilcox[Wilcoxon tests\nsigned -log10(p)\nBonferroni threshold]
  qcSig{miRNA sig in >= 2 clusters\nand let-7 non-sig?}
  tuneTopN[Tune analysis\nadjust topN / exclude artefacts\nrecompute stats]

  viz[Visualize & export\nheatmaps / surface plots\nHTML/PNGs]
  shiny[Shiny app\nload .RDS -> interactive plots]
  fin((End))

  %% Main spine
  start --> input --> validate --> reqOK
  reqOK -- Yes --> init --> ae --> denoise --> cluster --> freeze --> enoughK
  enoughK -- Yes --> targets --> logfc --> wilcox --> qcSig
  qcSig -- Yes --> viz --> shiny --> fin

  %% Short loops
  reqOK -- No --> fixInputs --> validate
  enoughK -- No --> tuneClust --> cluster
  qcSig -- No --> tuneTopN --> targets
```
## Getting Starteda

Follow these steps to set up the project environment and prepare the pipeline for use:

1. **Prerequisites – R Setup:** Install **R** (the pipeline was developed and tested on R 4.x) and optionally RStudio for a convenient development environment. Make sure you have internet access from R to install packages. If on Windows, also install Rtools (for compiling packages from source if needed).

2. **Clone the Repository:** Download or clone this repository to your local machine. For example, via command line:  
   ```bash
   git clone https://github.com/your-username/10X-Visium-Spatial-Transcriptomics-Pipeline-to-Analyse-miRNA-Expression-Patterns.git
   ```  
   Then set your R working directory to the project folder.

3. **Install Required Packages:** The pipeline relies on several R packages (from CRAN, Bioconductor, and GitHub). Below is a list of the main packages and how to install them:
   ```r
   # In an R session, install CRAN packages
   install.packages(c("Seurat", "SeuratData", "ggplot2", "patchwork", 
                      "tidyverse", "dplyr", "stringr", "readr", 
                      "glmGamPoi", "gridExtra", "pheatmap"))
   # Install Bioconductor packages
   install.packages("BiocManager")
   BiocManager::install(c("BayesSpace", "spacexr"))
   # Install SPATA2 (version 2.0.4) from GitHub
   install.packages("devtools")
   devtools::install_github("theMILOlab/SPATA2", ref = "v2.0.4")
   # Install the R Keras package and its TensorFlow backend (for autoencoder)
   install.packages("keras")
   keras::install_keras()  # This will download and install TensorFlow
   ```  
   **Note:** It’s important to use the package versions above (for example, Seurat v5.0.3 and SPATA2 v2.0.4) to ensure compatibility with the code. Newer versions may have breaking changes. The `keras::install_keras()` step will set up a Python backend for the autoencoder; you may skip it if you already have a working Keras/TensorFlow installation.

4. **Obtain a Visium Dataset:** You will need a 10x Genomics Visium spatial transcriptomics dataset to analyze. If you want to replicate the thesis analysis, you should download the specific example datasets (e.g. from GEO or the original publications) used in the project – these include **heart** and **brain** tissue sections with known tissue-specific miRNAs. (For instance, the thesis utilized data from GEO accessions like GSM5691527, GSM5691529, etc.) Alternatively, you can use **your own Visium dataset**. Make sure you have the **Space Ranger output** directory for your sample (which contains the `filtered_feature_bc_matrix` and `spatial` folder with tissue image and spot coordinates). Place your data in a known location or inside the project directory.  

5. **Configure data paths:** If using your own data or the downloaded example data, update the file paths in the analysis scripts accordingly. In the R Markdown scripts (especially the initialization script), set the path to your Visium data (the folder containing the `outs` or the matrix files). By default, the code may refer to specific file names (e.g., `GSMxxxx_sampleXYZ/outs`); you will need to replace those with your actual paths. Similarly, ensure you have the TargetScan miRNA target data file (a CSV of predicted targets) if you plan to run the miRNA analysis – the code expects it at `data/Targetscan_Data/...csv`. (If you don’t have this file, you can download target prediction data from [TargetScan](https://www.targetscan.org/) or use the one provided in the thesis.) Create a `data` directory in the project if it doesn’t exist, and put the necessary files there, or adjust the paths in the scripts to where your data resides.

Once the environment is set up and data is in place, you can proceed to run the analysis as described below.

## Usage

### Running the Analysis Pipeline

The analysis is organized into a series of R Markdown (`.Rmd`) scripts in the `R_Markdown_Scripts` directory. The main entry point is **`Master_Analysis_Script_V4.Rmd`**, which ties together the full pipeline (from data processing to generating outputs). To run the pipeline:

- **Interactive (RStudio):** Open the `Master_Analysis_Script_V4.Rmd` in RStudio. Ensure the working directory is set to the project folder. You might need to run the preliminary setup code in `Object_initialisation_V4.Rmd` first if you are starting from raw data. This will create the SPATA2 object and perform normalization/denoising for your dataset, saving intermediate results (e.g., as RDS files). After that, run (knit) the `Master_Analysis_Script_V4.Rmd` notebook. This will execute all sections of the pipeline in order. It will produce output files including interactive HTML reports (saved under `HTML_outputs/`) and plots illustrating the spatial expression patterns and cluster analyses.

- **Batch/Scripted:** Alternatively, you can run the pipeline non-interactively. For example, from an R console you could call:  
  ```r
  rmarkdown::render("R_Markdown_Scripts/Master_Analysis_Script_V4.Rmd")
  ```  
  Make sure to edit any file paths in the scripts to point to your data before running. The `Master_Analysis_Script_V4.Rmd` assumes that the SPATA2 object (with denoised data) is ready to be loaded from an RDS file (to save time on repeated runs). If you haven’t created this yet, run the `Object_initialisation_V4.Rmd` to generate it. After a successful run, check the `HTML_outputs` folder for the results – you can open the `.html` files in a web browser to interactively explore the heatmaps and tables of results.

The pipeline will output various results, such as:
- **Cluster tables:** Lists of top marker genes for each cluster and statistical comparisons between clusters.
- **miRNA target overlaps:** Tables showing which predicted miRNA targets are enriched in each cluster’s marker set or at cluster boundaries.
- **Heatmaps and plots:** Spatial heatmaps (in HTML reports) highlighting gene expression differences (e.g., each cluster vs the rest, or pairwise cluster comparisons with significance indicated). These help visualize potential miRNA effects at the boundaries.
- If everything is set up correctly, the R Markdown should complete without errors and produce a comprehensive report of the analysis.

### Launching the Shiny App

The repository includes a Shiny application (in the `Shiny_App_Script` folder, file `SPATA2_app.R`) for interactive exploration of the spatial transcriptomics data and results. To launch the app:

1. Ensure you have completed the pipeline analysis for at least one dataset and have the resulting SPATA2 object saved as an RDS file (or use the example RDS if provided). By default, the app looks for `.rds` files in a folder `data/RDS_Files_V3_Script/` (this can be changed by editing `rds_directory` at the top of `SPATA2_app.R`). Place your SPATA2 object file in that directory or update the path.

2. In an R session, make sure your working directory is the repository (or set it to the `Shiny_App_Script` directory). Then run:  
   ```r
   library(shiny)
   runApp("Shiny_App_Script")
   ```  
   This will start the Shiny app locally. (In RStudio, you can also click the **Run App** button with `SPATA2_app.R` open.)

3. Once the app is running, you’ll see a web interface open. In the app:
   - Use the dropdown to **select a SPATA2 RDS file** (if you placed one in the data directory, it should appear in the list). This will load the spatial transcriptomics data for that sample.
   - The app provides controls to choose a cluster of interest and a comparison cluster (or all other cells) to generate heatmaps. You can select which cluster’s markers to visualize and how to compare clusters.
   - You can also select **miRNAs or gene sets** if the app includes those options, to highlight where their target genes are expressed.
   - The output will show an interactive heatmap or spatial plot in the main panel, which you can hover for details or download.

4. Use the Shiny app to explore different clusters and boundaries. For example, you might select cluster A vs cluster B to see genes enriched in A relative to B, and observe if those genes include targets of a specific miRNA. The app is a convenient way to test various scenarios without rerunning the R Markdown each time.

**Note:** The Shiny app is for exploration and uses the processed data; any substantial changes to the analysis (e.g., re-running clustering with different parameters) should be done through the R scripts and then re-loaded into the app.

## Acknowledgments

This pipeline was developed as part of an **MSc Bioinformatics project at the University of Edinburgh (2023–2024)**. The work was inspired by the hypothesis that miRNAs act as regulators of gene expression at the interface of different cell types. We acknowledge the support of the University of Edinburgh and the project supervisor throughout the development of this tool. The repository and code are made available for educational and research purposes, and we hope they prove useful to others investigating spatial transcriptomics and miRNA regulation.
