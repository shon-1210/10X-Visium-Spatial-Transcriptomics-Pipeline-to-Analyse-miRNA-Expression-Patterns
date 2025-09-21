# Project Title: Using spatial transcriptomics to evaluate if microRNAs regulate expression profiles at cell-type boundaries 

## A 10X Visium Spatial Transcriptomics Pipeline to Analyse miRNA Expression Patterns

**Abstract:**  

MicroRNAs (miRNAs) regulate messenger RNA (mRNA) expression post-transcriptionally in mammals by non-slicing repression. When a miRNA forms the RNA-induced silencing complex (RISC) with an Argonaute (AGO) protein, it targets specific mRNAs which are known as its targets. Since their discovery in 1993 by Lee et al. in *C. elegans*, the functions of miRNAs have been widely debated, with various mechanisms proposed. Stark et al. found in *Drosophila* that certain miRNAs and their targets are not co-expressed but are instead found in separate cells, suggesting a regulatory decision at specific transition points. This project investigates whether miRNAs act as guardians of gene expression between neighbouring cell types by analysing the expression levels of their targets at these transition points.  

To explore this, we developed a pipeline to detect the effect of tissue-specific miRNAs in Spatial Transcriptomic (ST) datasets, where gene expression levels are mapped to tissue locations. Selected ST datasets were then analysed to find neighbouring cell types with high and low expression of miRNA targets, indicating potential miRNA activity. By visualising the spatial expression patterns of the targets at these transition points using ST analysis tools, evidence can be collected to support the hypothesis. miR-124 and miR-1, specifically expressed in Brain and Heart tissue respectively, were detected with high confidence and with neighbouring clusters of high and low expression obtained for both. However, the hypothesis could not be validated in the brain dataset that was analysed in detail, as the expression of the targets was low only in one region and expressed to a higher level in varying degrees for all other cellular regions. The results however were promising as an initial step and should be followed up on by testing more datasets.

## Key Features

- **End-to-end spatial analysis pipeline:** Process raw 10x Visium data (gene–barcode matrices and tissue images) into analysis-ready form using Seurat and SPATA2 frameworks.  
- **Denoising with autoencoders:** Improve data quality by applying an autoencoder to reduce technical noise in the expression matrix.  
- **Flexible clustering:** Identify spatial domains using Hartigan–Wong K-means and BayesSpace (spatial prior).  
- **miRNA target integration:** Use TargetScan-predicted miRNA–mRNA interactions to test whether differences at cluster interfaces are enriched for targets of tissue miRNAs.  
- **Spatial visualisation:** Generate heatmaps (targets vs non-targets), density plots, spatial trajectories, and cluster maps to inspect spatial patterns.  
- **Interactive Shiny app:** Explore the pipeline outputs interactively by selecting clustering/expression-matrix combinations and adjusting **topN** for target lists; reproduce the study heatmaps.

## How It Works

1. **Data Input & Setup:** Load 10x Visium outputs (filtered matrix and spatial image/coordinates) and initialise a SPATA2/Seurat object. Preprocess with SCTransform v2.  
2. **Denoising with Autoencoders:** Train an autoencoder and use the reconstruction to create a **denoised** expression layer for downstream analysis.  
3. **Clustering (K-means & BayesSpace):** Run BayesSpace and Hartigan–Wong K-means; set K-means **k** to the BayesSpace cluster count.  
4. **miRNA Target Analysis:** Compare targets vs non-targets per cluster and between neighbouring clusters to infer miRNA presence/absence.  
5. **Outputs:** Save RDS objects and dataset-specific HTML heatmaps for inspection; use the Shiny app for interactive exploration.

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
1. **Start** — Begin pipeline  
2. **10x Visium Data** — Load `filtered_feature_bc_matrix.h5` + `spatial/` folder  
3. **Validate Files** — Check paths, image, and gene–barcode matrix  
4. **Files Valid?** — Quality check decision point  
5. **Fix Paths** — Repair issues or re-run Space Ranger (if needed)

---

## 🔧 Section 2: Preprocessing

```
Initialize Object → Autoencoder Assessment → Denoise Matrix → PROCEED TO SECTION 3
```

**Steps:**
1. **Initialize Object** — Create SPATA2/Seurat object with SCTransform v2  
2. **Autoencoder Assessment** — Choose activation and bottleneck settings  
3. **Denoise Matrix** — Create the denoised expression layer

---

## 🎯 Section 3: Clustering + Save .RDS for downstream analysis

```
Run Clustering → Save Results → ≥3 Clusters?
                                     ↓
                                    NO → Tune Parameters → (back to Run Clustering)
                                     ↓
                                    YES → PROCEED TO SECTION 4
```

**Steps:**
1. **Run Clustering** — Run BayesSpace and K-means clustering techniques; set K-means **k** to match BayesSpace  
2. **Save Results** — Write the `.rds` object (with clustering and denoised data)  
3. **≥3 Clusters?** — If fewer than 3 distinct clusters of interest are found (or quality is poor), discard the dataset

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
1. **Load Targets** — Import TargetScan; choose **topN** (e.g., 100/200/300) and include **let-7** as a negative control  
2. **Expression Matrix** — Build the expression matrix from **log1p-transformed corrected counts**  
3. **Compute logFC** — Calculate logFC for cluster-vs-rest and pairwise neighbour-vs-neighbour comparisons (for each clustering method)  
4. **Statistical Tests** — Wilcoxon rank-sum tests comparing targets vs non-targets; **Bonferroni** threshold at **|log10 p| ≥ 4**  
5. **QC Pass?** — Require the tissue miRNA to be significant in **≥2 clusters** and **let-7** to be non-significant  
6. **Tune Analysis** — Adjust **topN**, revisit clustering to exclude artefacts, and recompute if needed

---

## 📊 Section 5: Visualisation & Output

```
Generate Plots → Shiny App → Final Outputs → PIPELINE COMPLETE!
```

**Steps:**
1. **Generate Plots** — Create **heatmaps** (log10-transformed p-values, targets vs non-targets), **density plots**, **spatial trajectories**, and cluster maps  
2. **Shiny App** — Reproduce the analysis heatmaps by selecting clustering/expression-matrix combinations and adjusting **topN**  
3. **Final Outputs** — Save pipeline `.rds` objects and dataset-specific HTML heatmaps (e.g., `Brain2_test_targets.html`, `Heart21d_test_targets.html`)

---

## 🔄 Pipeline Flow Summary

```
Section 1 → Section 2 → Section 3 → Section 4 → Section 5
   ↓           ↓           ↓           ↓           ↓
Data Input  Preprocess   Clustering   miRNA     Visualise
& Validate               + logFC      Analysis  & Output
```

## 📋 Quick Reference

**Input files required**
- `filtered_feature_bc_matrix.h5`  
- `spatial/` folder (image + spot coordinates)

**Key checkpoints**
- ✅ Files validate successfully  
- ✅ At least 3 distinct clusters  
- ✅ **logFC** from **log1p-transformed corrected counts**  
- ✅ Significance: **|log10 p| ≥ 4** (Bonferroni)  
- ✅ Tissue miRNA significant in **≥2 clusters**  
- ✅ **let-7** negative control is non-significant

**Final outputs**
- 📄 **`.rds` objects** — saved SPATA2 objects after clustering/denoising  
- 📈 **HTML heatmaps** — dataset-specific (e.g., `Brain1_test_targets.html`, `Brain2_test_targets.html`, `Heart21d_test_targets.html`)  
- 💻 **Shiny app** — interactive reproduction of the study heatmaps

**Technologies used**
- SPATA2 / Seurat for spatial data handling and preprocessing  
- BayesSpace and Hartigan–Wong K-means for clustering  
- TargetScan for miRNA target lists  
- Shiny for interactive visualisation

---

### Interactive Shiny App

**🌐 Live demo:** https://shonkuriangeorge.com/spata2shinyapp/  

The dissertation references an internal University of Edinburgh deployment; the link above is my public deployment for demonstration.  

**How to use:**
1. Select a **clustering/expression-matrix** combination  
2. Set **topN** for the miRNA target list and run the comparisons  
3. Inspect the resulting **heatmaps** (log10 p-values for targets vs non-targets) across clusters  
4. Explore cluster interfaces to assess potential miRNA regulation patterns

*For developers:* Source code is in `Shiny_App_Script/SPATA2_app.R`. Use pipeline-produced `.rds` objects to load your own processed datasets locally.

---

*This pipeline processes 10x Visium spatial transcriptomics data through five stages with built-in QC and parameter tuning at each step.*

## Acknowledgments

This pipeline was developed as part of an **MSc Bioinformatics project at the University of Edinburgh (2023–2024)**. The repository and code are provided for educational and research use. The complete thesis can be shared on request. Please contact shon.kurian.george@gmail.com if you wish to do so.

