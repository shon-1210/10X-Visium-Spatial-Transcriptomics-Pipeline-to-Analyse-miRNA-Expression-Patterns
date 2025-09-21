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

1. **Data Input & Setup:** The pipeline starts by loading a spatial transcriptomics dataset (10x Visium output, including the count matrix and spatial coordinates/image) and initializing it as a SPATA2/Seurat object. Basic preprocessing like quality control and normalization (using Seurat's SCTransform v2) is applied to prepare the data for analysis.  
2. **Denoising with Autoencoders:** An autoencoder neural network is then trained on the expression data to learn a compressed representation. The network's output is used to reconstruct a denoised expression matrix, reducing technical noise and emphasizing true biological signal.  
3. **Clustering (K-means & BayesSpace):** The pipeline performs unsupervised clustering on the spatial spots in two ways: (a) **K-means clustering** on the gene expression features to group spots by transcriptional similarity, and (b) **BayesSpace clustering**, which leverages a Bayesian model and spatial information to group neighboring spots into the same cluster. These approaches help identify distinct cell-type regions or spatial domains within the tissue.  
4. **miRNA Target Analysis:** For each identified cluster and boundary region, the pipeline analyzes gene expression differences and cross-references them with known miRNA target data. Using TargetScan predictions, it checks whether the genes up- or down-regulated at cluster edges are enriched for targets of certain miRNAs (especially tissue-specific miRNAs). This helps evaluate the hypothesis that those miRNAs might be regulating genes at the interface of different cell types.  
5. **Visualization & Outputs:** The results are compiled into easily interpretable outputs. The pipeline generates tables of cluster-specific marker genes and miRNA target overlaps, as well as **spatial heatmaps** that illustrate where particular gene expression or miRNA activity is high or low across the tissue slice. For example, it can produce a heatmap showing a cluster's top genes compared to all other clusters, or pairwise cluster comparison plots highlighting boundary effects. These visualizations are output as HTML reports and/or image files for further inspection.  
6. **Interactive Exploration:** Finally, the included Shiny app allows users to interactively explore the processed data. Users can load a processed SPATA2 object (as an `.rds` file), then dynamically generate plots and compare clusters via a web interface – without needing to rerun code. This is especially useful for examining specific genes or miRNAs of interest in the spatial context, beyond the static results.  

## Pipeline Flowchart

Below is a high-level flowchart of the pipeline:

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

## 🎯 Section 3: Clustering + Save .RDS file for downstream analysis

```
Run Clustering → Save Results → ≥3 Clusters?
                                     ↓
                                    NO → Tune Parameters → (back to Run Clustering)
                                     ↓
                                    YES → PROCEED TO SECTION 4
```

**Steps:**
1. **Run Clustering** - Execute both BayesSpace & K-means clustering methods, set K-means number of clusters to match BayesSpace
2. **Save Results** - Export .RDS file with complete object
3. **≥3 Clusters?** - Check if a sufficient number of clusters are identified (accurately too), if not, discard the dataset


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
2. **Create Expression Matrix** - Build Expression Matrix from log1p() transformed Corrected Count Matrix
3. **Compute logFC** -  Calculate logFC as required for individual cluster vs rest comparisons and cluster vs cluster pairwise comparisons, done iteratively for each clustering method
4. **Statistical Tests** - Wilcoxon rank-sum statistical test performed using logFC values of targets and non-targets of the miRNA being analysed; Bonferroni correction to determine p-value significance threshold
5. **QC Pass?** - Check: miRNA significant in ≥2 clusters AND let-7 non-significant
6. **Tune Analysis** - Adjust topN, exclude artifacts, recompute stats (if needed)

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

### Interactive Shiny App

**🌐 Try the live app: [https://shonkuriangeorge.com/spata2shinyapp/](https://shonkuriangeorge.com/spata2shinyapp/)**

The pipeline includes an interactive Shiny web application for exploring spatial transcriptomics data and miRNA analysis results. You can access the app directly through the website above, which includes example datasets ready for exploration.

**How to use the app:**

1. **Load Dataset** - Use the dropdown menu to select a pre-loaded SPATA2 dataset (brain, heart, or other tissue samples)

2. **Explore Clusters** - Choose clusters of interest to compare:
   - Select a primary cluster and comparison cluster (or "all other cells")
   - View cluster-specific marker genes and their spatial distribution

3. **miRNA Analysis** - Investigate miRNA target enrichment:
   - Select specific miRNAs to examine their predicted targets
   - Compare target gene expression between neighboring clusters
   - Visualize spatial patterns of miRNA activity

4. **Interactive Visualization** - Generate and explore plots:
   - Interactive heatmaps showing gene expression patterns
   - Spatial plots overlaying expression data on tissue images
   - Hover for detailed information and download results

5. **Boundary Analysis** - Examine cluster interfaces:
   - Compare gene expression at cluster boundaries
   - Identify potential miRNA regulation zones
   - Test the hypothesis of miRNAs as spatial expression guardians

**For developers:** The complete Shiny app source code is available in the `Shiny_App_Script/` folder (file: `SPATA2_app.R`) if you want to run the app locally or modify it for your own datasets.

**Note:** The web app uses pre-processed example datasets. To analyze your own data, process it through the pipeline first, then either use the local version of the app or contact us about uploading custom datasets.

---

*This pipeline processes 10X Visium spatial transcriptomics data through five main stages with built-in quality control and parameter tuning at each step.*

## Acknowledgments

This pipeline was developed as part of an **MSc Bioinformatics project at the University of Edinburgh (2023–2024)**. The work was inspired by the hypothesis that miRNAs act as regulators of gene expression at the interface of different cell types. We acknowledge the support of the University of Edinburgh and the project supervisor throughout the development of this tool. The repository and code are made available for educational and research purposes, and we hope they prove useful to others investigating spatial transcriptomics and miRNA regulation.
