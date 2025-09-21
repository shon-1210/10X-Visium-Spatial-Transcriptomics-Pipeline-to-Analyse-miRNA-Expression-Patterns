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

# 10X Visium miRNA Spatial Transcriptomics Pipeline
*Mobile-Friendly Version*

## 📁 Section 1: Data Input & Validation

```
START → 10x Visium Data → Validate Files → Files Valid?
                                              ↓
                                             NO → Fix Paths → (back to Validate)
                                              ↓
                                             YES → PROCEED TO SECTION 2
```

**Steps:**
1. **Start** - Begin pipeline
2. **10x Visium Data** - Load filtered_feature_bc_matrix.h5 + spatial folder
3. **Validate Files** - Check paths, image, and gene-barcode matrix
4. **Files Valid?** - Quality check decision point
5. **Fix Paths** - Repair issues or re-run Space Ranger (if needed)

---

## 🔧 Section 2: Preprocessing

```
Initialize Object → Autoencoder Assessment → Denoise Matrix → PROCEED TO SECTION 3
```

**Steps:**
1. **Initialize Object** - Create SPATA2/Seurat object with SCTransform v2
2. **Autoencoder Assessment** - Select optimal activation and bottleneck parameters
3. **Denoise Matrix** - Create denoised expression layer for analysis

---

## 🎯 Section 3: Spatial Clustering

```
Run Clustering → Save Results → ≥3 Clusters?
                                     ↓
                                    NO → Tune Parameters → (back to Run Clustering)
                                     ↓
                                    YES → PROCEED TO SECTION 4
```

**Steps:**
1. **Run Clustering** - Execute BayesSpace + K-means with spatial prior
2. **Save Results** - Export .RDS file with complete object
3. **≥3 Clusters?** - Check if sufficient clusters identified
4. **Tune Parameters** - Adjust BayesSpace q parameters, re-run K-means (if needed)

---

## 🧬 Section 4: miRNA Target Analysis

```
Load Targets → Compute logFC → Statistical Tests → QC Pass?
                                                      ↓
                                                     NO → Tune Analysis → (back to Load Targets)
                                                      ↓
                                                     YES → PROCEED TO SECTION 5
```

**Steps:**
1. **Load Targets** - Import TargetScan database, choose topN + let-7 control
2. **Compute logFC** - Calculate cluster vs rest, pairwise comparisons
3. **Statistical Tests** - Wilcoxon rank-sum with Bonferroni correction
4. **QC Pass?** - Check: miRNA significant in ≥2 clusters AND let-7 non-significant
5. **Tune Analysis** - Adjust topN, exclude artifacts, recompute stats (if needed)

---

## 📊 Section 5: Visualization & Output

```
Generate Plots → Shiny App → Final Outputs → PIPELINE COMPLETE! 🎉
```

**Steps:**
1. **Generate Plots** - Create heatmaps and surface plots
2. **Shiny App** - Build interactive visualization interface
3. **Final Outputs** - Export analysis.RDS, plots.html, Shiny interface

---

## 🔄 Pipeline Flow Summary

```
Section 1 → Section 2 → Section 3 → Section 4 → Section 5
   ↓           ↓           ↓           ↓           ↓
Data Input  Preprocess  Clustering   miRNA     Visualize
& Validate              Spatial    Analysis   & Output
```

## 📋 Quick Reference

**Input Files Required:**
- `filtered_feature_bc_matrix.h5`
- `spatial/` folder with images and coordinates

**Key Quality Checkpoints:**
- ✅ File validation
- ✅ Minimum 3 clusters
- ✅ miRNA significance in ≥2 clusters
- ✅ let-7 control non-significant

**Final Outputs:**
- 📄 `analysis.RDS` - Complete results object
- 📈 `plots.html` - Static visualizations  
- 💻 Shiny app - Interactive interface

**Technologies Used:**
- SPATA2/Seurat for spatial analysis
- BayesSpace for spatial clustering
- TargetScan for miRNA targets
- Shiny for interactive visualization

---

*This pipeline processes 10X Visium spatial transcriptomics data through five main stages with built-in quality control and parameter tuning at each step.*

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
