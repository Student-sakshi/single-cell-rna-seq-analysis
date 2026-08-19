
# Sakshi Parate - 7073022

suppressPackageStartupMessages({
  library(dplyr)
  library(spatstat.core)
  library(Seurat)
  library(patchwork)
  library(DoubletFinder)
  library(SingleR)
  library(enrichR)
  library(SingleCellExperiment)
  library(SeuratWrappers)
  library(tidyverse)
})

# Load all four RDS files
sample1 <- readRDS("/Users/sam/Documents/UdS/Sem3/SCB/Assignments/Project1/dataset/GSM4557329_GSM4557329_556_cell.counts.matrices.rds")
sample2 <- readRDS("/Users/sam/Documents/UdS/Sem3/SCB/Assignments/Project1/dataset/GSM4557330_GSM4557330_557_cell.counts.matrices.rds")
sample3 <- readRDS("/Users/sam/Documents/UdS/Sem3/SCB/Assignments/Project1/dataset/GSM4557331_GSM4557331_558_cell.counts.matrices.rds")
sample4 <- readRDS("/Users/sam/Documents/UdS/Sem3/SCB/Assignments/Project1/dataset/GSM4557337_GSM4557337_HIP043_cell.counts.matrices.rds")

# Create a Seurat object for each sample
seurat1 <- CreateSeuratObject(counts = sample1, project = "covid_556", min.cells = 3, min.features = 200)
seurat2 <- CreateSeuratObject(counts = sample2, project = "covid_557", min.cells = 3, min.features = 200)
seurat3 <- CreateSeuratObject(counts = sample3, project = "covid_558", min.cells = 3, min.features = 200)
seurat4 <- CreateSeuratObject(counts = sample4, project = "HIP043",    min.cells = 3, min.features = 200)

# Label each sample with the corresponding metadata
sample_metadata <- data.frame(
  Sample = c("covid_556", "covid_557", "covid_558", "HIP043"),
  Donor = c("C2", "C3", "C4", "H4"),
  Replicate = c("T1", "T1", "T1", "T1"),
  Sex = c("F", "M", "M", "F")
)
sample_metadata

# Add metadata to each Seurat object
seurat1@meta.data$Donor <- "C2"
seurat1@meta.data$Replicate <- "T1"
seurat1@meta.data$Sex <- "F"

seurat2@meta.data$Donor <- "C3"
seurat2@meta.data$Replicate <- "T1"
seurat2@meta.data$Sex <- "M"

seurat3@meta.data$Donor <- "C4"
seurat3@meta.data$Replicate <- "T1"
seurat3@meta.data$Sex <- "M"

seurat4@meta.data$Donor <- "H4"
seurat4@meta.data$Replicate <- "T1"
seurat4@meta.data$Sex <- "F"

# number of cells
ncol(seurat1)
ncol(seurat2)
ncol(seurat3)
ncol(seurat4)

#number of genes
nrow(seurat1)
nrow(seurat2)
nrow(seurat3)
nrow(seurat4)

#metadata information
colnames(seurat1@meta.data)
head(seurat1@meta.data)

saveRDS(data, file ="scb.R")

# ---- 4: Preprocessing ----

# ---- 4.1: Preprocessing ----

# Add percent mitochondrial metric
seurat1[["percent.mt"]] <- PercentageFeatureSet(seurat1, pattern = "^MT-")
seurat2[["percent.mt"]] <- PercentageFeatureSet(seurat2, pattern = "^MT-")
seurat3[["percent.mt"]] <- PercentageFeatureSet(seurat3, pattern = "^MT-")
seurat4[["percent.mt"]] <- PercentageFeatureSet(seurat4, pattern = "^MT-")

# Apply filtering
seurat1 <- subset(seurat1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
seurat2 <- subset(seurat2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
seurat3 <- subset(seurat3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
seurat4 <- subset(seurat4, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Doublet Finder
# Combine samples into a list
seurat_list <- list(
  seurat1 = seurat1,
  seurat2 = seurat2,
  seurat3 = seurat3,
  seurat4 = seurat4
)

df_results <- list()
expected_rate <- 0.06   # 6% expected doublet rate

for (s in names(seurat_list)) {
  
  message("---- Processing ", s, " ----")
  obj <- seurat_list[[s]]
  
  # 1. Preprocessing needed for DoubletFinder
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = 20)
  obj <- FindNeighbors(obj, dims = 1:10)
  obj <- FindClusters(obj, resolution = 0.5)
  
  # 2. Parameter sweep to find optimal pK
  sweep.res <- paramSweep(obj, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res)
  pk_stats <- find.pK(sweep.stats)
  
  best.pK <- as.numeric(as.character(
    pk_stats$pK[which.max(pk_stats$BCmetric)]
  ))
  
  message("Best pK for ", s, ": ", best.pK)
  
  # 3. Estimate expected doublets
  nCells <- ncol(obj)
  nExp <- round(nCells * expected_rate)
  
  # Adjust for homotypic doublets
  homotypic.prop <- modelHomotypic(obj$seurat_clusters)
  nExp.adj <- round(nExp * (1 - homotypic.prop))
  
  # 4. Run DoubletFinder
  obj <- doubletFinder(
    obj,
    PCs = 1:10,
    pN = 0.25,
    pK = best.pK,
    nExp = nExp.adj
  )
  
  # 5. Identify DoubletFinder classification column
  df_col <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
  df_col <- tail(df_col, 1)   # take the latest generated column
  
  message("DF column: ", df_col)
  
  # 6. Remove doublets 
  obj_clean <- obj[, obj@meta.data[[df_col]] == "Singlet"]
  message("Remaining singlet cells: ", ncol(obj_clean))
  
  # Store cleaned object
  seurat_list[[s]] <- obj_clean
  
  # Save summary info
  df_results[[s]] <- list(
    pK = best.pK,
    nExp = nExp,
    nExp_adjusted = nExp.adj,
    DF_column = df_col,
    cells_after = ncol(obj_clean)
  )
}

# Extract cleaned objects
seurat1 <- seurat_list$seurat1
seurat2 <- seurat_list$seurat2
seurat3 <- seurat_list$seurat3
seurat4 <- seurat_list$seurat4

df_results

# Normalization and Feature selection
# List of cleaned Seurat objects
seurat_list <- list(
  seurat1 = seurat1,
  seurat2 = seurat2,
  seurat3 = seurat3,
  seurat4 = seurat4
)

# Apply normalization and feature selection to each
for (s in names(seurat_list)) {
  obj <- seurat_list[[s]]
  
  # Normalization (default: LogNormalize)
  obj <- NormalizeData(obj)
  
  # Find variable features
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  
  # Save back
  seurat_list[[s]] <- obj
}

# Extract updated objects
seurat1 <- seurat_list$seurat1
seurat2 <- seurat_list$seurat2
seurat3 <- seurat_list$seurat3
seurat4 <- seurat_list$seurat4

# ---- 4.2A: MERGING WITHOUT BATCH CORRECTION ----

merged_noBC <- merge(
  seurat1,
  y = list(seurat2, seurat3, seurat4),
  add.cell.ids = c("556", "557", "558", "HIP043"),
  project = "PBMC_noBatchCorrection"
)

# Normalize + PCA + UMAP for visualization
merged_noBC <- NormalizeData(merged_noBC)
merged_noBC <- FindVariableFeatures(merged_noBC)
merged_noBC <- ScaleData(merged_noBC)
merged_noBC <- RunPCA(merged_noBC, npcs = 30)
merged_noBC <- RunUMAP(merged_noBC, dims = 1:20)

# Plot UMAP colored by sample to see batch effect
DimPlot(merged_noBC, group.by = "orig.ident", label = TRUE) +
  ggtitle("A: UMAP WITHOUT Batch Correction")

# --- 4.2B: SEURAT DATA INTEGRATION (WITH BATCH CORRECTION) ----

# Put cleaned objects into list
seurat_list <- list(seurat1, seurat2, seurat3, seurat4)

# 1. Normalize and find variable features 
for (i in 1:length(seurat_list)) {
  seurat_list[[i]] <- NormalizeData(seurat_list[[i]])
  seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]], nfeatures = 2000)
}

# 2. Select integration features
features <- SelectIntegrationFeatures(seurat_list)

# 3. Find integration anchors
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  anchor.features = features,
  dims = 1:30
)

# 4. Integrate
integrated <- IntegrateData(
  anchorset = anchors,
  dims = 1:30
)

# 5. Scaling + PCA + UMAP
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, npcs = 30)
integrated <- RunUMAP(integrated, dims = 1:20)

# 6. Clustering
integrated <- FindNeighbors(integrated, dims = 1:20)
integrated <- FindClusters(integrated, resolution = 0.5)

# 7. Plot
DimPlot(integrated, group.by = "orig.ident") +
  ggtitle("UMAP WITH Batch Correction (LogNormalize Integration)")
  
# --- 5: Dimensionality Reduction ---

# --- 5.1: Dimensionality Reduction ---

# PCA
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)

# Visualize PCA
DimPlot(integrated, reduction = "pca", group.by = "orig.ident") +
  ggtitle("PCA Colored by Sample")

# Elbow plot to select number of PCs
ElbowPlot(integrated, ndims = 50) +
  ggtitle("Elbow Plot for Selecting Number of PCs")

# Choose PCs based on Elbow 
pcs_to_use <- 1:20

# UMAP 
integrated <- RunUMAP(integrated, dims = pcs_to_use)

# UMAP plot
DimPlot(integrated, reduction = "umap", group.by = "orig.ident") +
  ggtitle("UMAP after PCA")

# --- 5.2: Clustering ---

# Rerun neighbors with selected PCs
integrated <- FindNeighbors(integrated, dims = pcs_to_use)

# Louvain
integrated <- FindClusters(integrated, resolution = 0.5, algorithm = 1)  
Idents(integrated) <- "seurat_clusters"
DimPlot(integrated, reduction = "umap", label = TRUE) +
  ggtitle("Louvain Clustering")

# Leiden
install.packages("leiden")
library("igraph")
library("leiden")

integrated <- FindNeighbors(integrated, dims = 1:20)
names(integrated@graphs)
g <- graph_from_adjacency_matrix(
  integrated@graphs$integrated_snn,
  mode = "undirected",
  diag = FALSE,
  weighted = TRUE
)
V(g)$name <- as.character(seq_len(vcount(g)))

# Run Leiden
leiden_clusters <- leiden(g, resolution_parameter = 0.5)
integrated$leiden <- as.factor(leiden_clusters)
Idents(integrated) <- "leiden"

# Plot
DimPlot(integrated, reduction = "umap", label = TRUE) +
  ggtitle("Leiden Clustering")

saveRDS(data, file ="scb.R")

# --- 6: Automatic Cell Type Annotation ---

# --- 6.1: Automatic Annotation ---

BiocManager::install("SingleR-inc/celldex")
library(celldex)

integrated[["integrated"]]
int_data <- GetAssayData(integrated, assay = "integrated", slot = "data")
genes <- rownames(int_data)
rna_counts <- round(exp(int_data) - 1)
rna_logdata <- int_data
sce <- SingleCellExperiment(
  assays = list(
    counts = rna_counts,
    logcounts = rna_logdata
  ),
  colData = integrated@meta.data
)
ref <- celldex::HumanPrimaryCellAtlasData()
assayNames(ref)
ref_data <- assay(ref, "logcounts")
ref_labels <- ref$label.main
test_data <- logcounts(sce)
getAnywhere(SingleR)
pred <- SingleR(
  sc_data  = test_data,
  ref_data = ref_data,
  types    = ref_labels,
  fine.tune = FALSE      # <- makes it fast         
)
integrated$SingleR_label <- pred$labels
DimPlot(integrated, reduction = "umap", group.by = "SingleR_label", label = TRUE)

# --- 6.2: Manual Annotation ---

library(dplyr)
library(Matrix)

# 1) Marker list 
marker.list <- list(
  HSC    = c("CD34", "CD38", "SCA1", "KIT"),
  LMPP   = c("CD38", "CD52", "CSF3R", "CA1", "KIT", "CD34", "FLK2"),
  CLP    = c("IL7R"),
  GMP    = c("ELANE"),
  CMP    = c("IL3", "CSF2", "CSF1"),   # GM-CSF -> CSF2, M-CSF -> CSF1 (use gene symbols)
  
  B      = c("CD19", "CD20", "CD38"),
  PreB   = c("CD19", "CD34"),
  Plasma = c("SDC1", "IGHA1", "IGLC1", "MZB1", "JCHAIN"),
  
  T      = c("CD3D"),
  CD8    = c("CD3D", "CD3E", "CD8A", "CD8B"),
  CD4    = c("CD3D", "CD3E", "CD4"),
  
  NK     = c("FCGR3A", "NCAM1", "NKG7", "KLRB1"),
  
  Ery    = c("GATA1", "HBB", "HBA1", "HBA2"),
  
  pDC    = c("IRF8", "IRF4", "IRF7"),
  cDC    = c("CD1C", "CD207", "ITGAM", "NOTCH2", "SIRPA"),
  
  Mono14 = c("CD14", "CCL3", "CCL4", "IL1B"),
  Mono16 = c("FCGR3A", "CD68", "S100A12"),
  
  Basophils = c("GATA2")
)

# 2) Prepare case insensitive marker list and dataset gene names
marker.list_lower <- lapply(marker.list, function(x) tolower(x))

# Use the RNA assay (full feature set) rather than the integrated assay (which may only contain 2000 features)
expr_mat <- GetAssayData(integrated, assay = "RNA", slot = "data") # log-normalised expression
gene_names <- rownames(expr_mat)
gene_names_lower <- tolower(gene_names)

# Check how many of each marker appear in the RNA assay
marker_hits_per_type <- sapply(marker.list_lower, function(m) sum(m %in% gene_names_lower))
print("marker hits per marker-set:")
print(marker_hits_per_type)

# 3) Prepare clusters and storage
Idents(integrated) <- integrated$seurat_clusters
cluster_ids <- levels(Idents(integrated))   # character vector "0"..."19"
cluster_to_celltype <- setNames(rep(NA_character_, length(cluster_ids)), cluster_ids)

# 4) Score clusters
# compute per cluster average using RNA assay
for (cl in cluster_ids) {
  cells <- WhichCells(integrated, idents = cl)
  if (length(cells) == 0) { next }
  # rowMeans on sparse matrix is fine
  avg <- Matrix::rowMeans(expr_mat[, cells, drop = FALSE])
  names(avg) <- gene_names  # ensure names are real-case
  
  # For each marker set compute mean(avg_expr of hits) 
  scores <- sapply(marker.list_lower, function(marker_genes_lower) {
    hits_idx <- which(gene_names_lower %in% marker_genes_lower)
    if (length(hits_idx) == 0) return(0)
    mean(avg[hits_idx], na.rm = TRUE)
  })
  
  # choose top scoring marker set
  best <- names(scores)[which.max(scores)]
  cluster_to_celltype[cl] <- best
}

print("cluster -> chosen marker set (abbrev):")
print(cluster_to_celltype)

# 5) Pretty mapping from abbrev to label
celltype.pretty <- c(
  HSC = "HSC (HSC)",
  LMPP = "LMPP (LMPP)",
  CLP = "CLP (CLP)",
  GMP = "GMP (GMP)",
  CMP = "CMP (CMP)",
  B = "B cell (B)",
  PreB = "Pre-B cell (PreB)",
  Plasma = "Plasma cell (Plasma)",
  T = "T cell (T)",
  CD8 = "CD8 T cell (CD8)",
  CD4 = "CD4 T cell (CD4)",
  NK = "NK cell (NK)",
  Ery = "Erythroblast (Ery)",
  pDC = "Plasmacytoid DC (pDC)",
  cDC = "Conventional DC (cDC)",
  Mono14 = "Monocyte CD14 (Mono14)",
  Mono16 = "Monocyte CD16 (Mono16)",
  Basophils = "Basophil"
)

# 6) Build final cluster 
cluster_pretty <- sapply(cluster_to_celltype, function(abbrev) {
  if (is.na(abbrev)) return("Unknown") #if cluster_to_celltype has unexpected abbreviations, fall back to "Unknown"
  lab <- celltype.pretty[[abbrev]]
  if (is.null(lab)) {
    # fallback to show the raw abbrev capitalized
    return(paste0(toupper(substring(abbrev,1,1)), substring(abbrev,2), " (", abbrev, ")"))
  } else {
    return(lab)
  }
}, USE.NAMES = FALSE)

manual_labels_map <- setNames(cluster_pretty, cluster_ids)  # named by cluster ids

# 7) Assign per cell labels to index by cluster names)
integrated$manual_celltype <- unname(manual_labels_map[ as.character(Idents(integrated)) ])

# sanity check
message("Manual annotation counts:")
print(table(integrated$manual_celltype, useNA = "ifany"))

# 8) Plot manual annotation UMAP
DimPlot(integrated, reduction = "umap", group.by = "manual_celltype", label = TRUE, repel = TRUE) +
  ggtitle("Manual Cell Type Annotation (Table 2 markers)")

# 9) Compare Manual vs Automatic
if ("SingleR_label" %in% colnames(integrated@meta.data)) {
  p1 <- DimPlot(integrated, group.by = "SingleR_label", label = TRUE) + ggtitle("Automatic")
  p2 <- DimPlot(integrated, group.by = "manual_celltype", label = TRUE) + ggtitle("Manual")
  print(p1 + p2)
}

# 10) Violin and UMAP for 3 genes
markers_to_plot <- c("CD3D","MS4A1","NKG7")
plots <- VlnPlot(integrated, features = markers_to_plot, group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots, ncol = 3) + plot_annotation(title = "Expression of 3 marker genes")
FeaturePlot(integrated, features = markers_to_plot, reduction = "umap")

# 11) Merge clusters into the cell type groups 
Idents(integrated) <- integrated$manual_celltype
integrated$celltype_merged <- Idents(integrated)
message("Finished manual annotation and merging.")

# --- 6.3: Cell Type Proportions ---

library(dplyr)
library(ggplot2)

Idents(integrated) <- integrated$manual_celltype # ensure manual annotation is active identity
integrated$sample <- integrated$orig.ident   # if sample names stored in 'orig.ident'

# 1. Compute proportions
prop_df <- integrated@meta.data %>%
  group_by(sample, manual_celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count))
print(prop_df)

# 2. Plot as barplot (stacked)
ggplot(prop_df, aes(x = sample, y = proportion, fill = manual_celltype)) +
  geom_bar(stat = "identity") +
  ylab("Proportion") +
  xlab("Sample") +
  ggtitle("Cell-Type Proportions per Sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- 7: Differential Expression Analysis ---

# --- 7.1: Differential Expression Analysis on cell-types ---

raw_counts <- GetAssayData(integrated, assay = "RNA", slot = "counts")

# find duplicated genes
dups <- duplicated(rownames(raw_counts))
message("Duplicate genes removed: ", sum(dups))

# keep only unique genes
raw_counts_clean <- raw_counts[!dups, ]

# rebuild RNA assay
newRNA <- CreateAssayObject(counts = raw_counts_clean)

# replace RNA assay
integrated[["RNA"]] <- newRNA

DefaultAssay(integrated) <- "RNA"

# Recreate normalized data (log-normalized)
integrated <- NormalizeData(integrated)
integrated <- FindVariableFeatures(integrated)

message("RNA assay cleaned and normalized successfully.")

Idents(integrated) <- integrated$manual_celltype
levels(Idents(integrated))

de_b_vs_t <- FindMarkers(
  integrated,
  ident.1 = "B cell (B)",
  ident.2 = "T cell (T)",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"
)
de_b_vs_t$gene <- rownames(de_b_vs_t)

write.csv(de_b_vs_t, "DE_B_vs_T.csv", row.names = FALSE)
message("B vs T DE complete.")

de_t_vs_mono <- FindMarkers(
  integrated,
  ident.1 = "T cell (T)",
  ident.2 = "Monocyte CD16 (Mono16)",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"
)
de_t_vs_mono$gene <- rownames(de_t_vs_mono)

write.csv(de_t_vs_mono, "DE_T_vs_Mono16.csv", row.names = FALSE)
message("T vs Mono16 DE complete.")

library(ggplot2)

# B vs T
volcano_b_t <- ggplot(de_b_vs_t, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-0.5, 0.5), color = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed") +
  ggtitle("Volcano Plot: B cells vs T cells") +
  xlab("log2 Fold Change") +
  ylab("-log10 adjusted p-value") +
  theme_bw()

print(volcano_b_t)

# T vs Monocyte
volcano_t_mono <- ggplot(de_t_vs_mono, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-0.5, 0.5), color = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed") +
  ggtitle("Volcano Plot: T cells vs Monocyte (CD16)") +
  xlab("log2 Fold Change") +
  ylab("-log10 adjusted p-value") +
  theme_bw()

print(volcano_t_mono)

# --- 7.1.1 Memory Formation on Cells Healthy vs Covid ---

# Create a condition label
integrated$condition <- ifelse(integrated$orig.ident == "HIP043", "Healthy", "COVID")
Idents(integrated) <- integrated$condition
levels(Idents(integrated))

# Run DE: COVID vs Healthy
de_covid_vs_healthy <- FindMarkers(
  integrated,
  ident.1 = "COVID",
  ident.2 = "Healthy",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"
)
de_covid_vs_healthy$gene <- rownames(de_covid_vs_healthy)

# Save results
write.csv(de_covid_vs_healthy, "DE_COVID_vs_Healthy.csv", row.names = FALSE)

# Volcano plot
library(ggplot2)

volcano_covid_healthy <- ggplot(de_covid_vs_healthy, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-0.5, 0.5), color = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed") +
  ggtitle("Volcano Plot: COVID vs Healthy") +
  xlab("log2 Fold Change") +
  ylab("-log10 adjusted p-value") +
  theme_bw()

print(volcano_covid_healthy)

top_up <- head(de_covid_vs_healthy[order(-de_covid_vs_healthy$avg_log2FC), ], 10)
top_down <- head(de_covid_vs_healthy[order(de_covid_vs_healthy$avg_log2FC), ], 10)
top_up
top_down
head(de_covid_vs_healthy[order(de_covid_vs_healthy$p_val_adj), ], 10)

# --- 7.2 Plot Differentially expressed genes ---

library(dplyr)
library(ggplot2)
library(tidyr)

# force DE result into a plain data.frame (removes Rle/list/factor issues)
clean_DE <- function(df) {
  # if matrix or other, coerce to data.frame first
  df <- as.data.frame(df)
  # move rownames into column 'gene'
  df$gene <- rownames(df)
  rownames(df) <- NULL
  # convert columns to atomic vectors
  for (col in colnames(df)) {
    # Rle -> vector
    if (inherits(df[[col]], "Rle")) df[[col]] <- as.vector(df[[col]])
    # lists -> unlist
    if (is.list(df[[col]])) df[[col]] <- unlist(df[[col]])
    # factor -> character
    if (is.factor(df[[col]])) df[[col]] <- as.character(df[[col]])
    # attempt to convert character numeric -> numeric
    if (is.character(df[[col]])) {
      suppressWarnings(num <- as.numeric(df[[col]]))
      if (!all(is.na(num))) df[[col]] <- num
    }
  }
  # final safe coercion
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  return(df)
}

# Clean your DE tables 
de_b_vs_t <- clean_DE(de_b_vs_t)
de_t_vs_mono <- clean_DE(de_t_vs_mono)

# Ensure p_val_adj column exists (Seurat uses p_val_adj); if not, try p_val
if (!"p_val_adj" %in% colnames(de_b_vs_t) & "p_val" %in% colnames(de_b_vs_t)) {
  de_b_vs_t$p_val_adj <- de_b_vs_t$p_val
}
if (!"p_val_adj" %in% colnames(de_t_vs_mono) & "p_val" %in% colnames(de_t_vs_mono)) {
  de_t_vs_mono$p_val_adj <- de_t_vs_mono$p_val
}

# Add derived columns and identify which cell type has higher expression
# For B vs T: positive avg_log2FC => higher in ident.1 (B), negative => higher in ident.2 (T)
de_b_vs_t <- de_b_vs_t %>%
  mutate(
    avg_log2FC = as.numeric(avg_log2FC), 
    p_val_adj = as.numeric(p_val_adj),
    neglog10_padj = -log10(pmax(p_val_adj, 1e-300)),
    enriched_in = ifelse(avg_log2FC > 0, "B cell (B)", "T cell (T)")
  )

# For T vs Mono16: positive => higher in ident.1 (T), negative => higher in ident.2 (Mono16)
de_t_vs_mono <- de_t_vs_mono %>%
  mutate(
    avg_log2FC = as.numeric(avg_log2FC),
    p_val_adj = as.numeric(p_val_adj),
    neglog10_padj = -log10(pmax(p_val_adj, 1e-300)),
    enriched_in = ifelse(avg_log2FC > 0, "T cell (T)", "Monocyte CD16 (Mono16)")
  )

# Select top 5 by adjusted p-value for each comparison
top5_b_t <- de_b_vs_t %>% arrange(p_val_adj) %>% slice_head(n = 5) %>% mutate(comparison = "B_vs_T")
top5_t_mono <- de_t_vs_mono %>% arrange(p_val_adj) %>% slice_head(n = 5) %>% mutate(comparison = "T_vs_Mono16")

# Combine
top5_all <- bind_rows(top5_b_t, top5_t_mono)

# For plotting, keep gene factor ordering per comparison so most significant are at top
top5_all <- top5_all %>%
  group_by(comparison) %>%
  mutate(gene = factor(gene, levels = rev(gene))) %>%
  ungroup()

# Combined dotplot: x = enriched cell type, y = gene, size = significance, color = fold-change
p_dot <- ggplot(top5_all, 
                aes(x = enriched_in, y = gene, size = neglog10_padj, color = avg_log2FC)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "log2FC") +
  scale_size_continuous(name = "-log10(adj p)") +
  facet_wrap(~ comparison, scales = "free_y", ncol = 1) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  ) +
  xlab("Cell type where the gene is enriched") +
  ylab("Top DE genes (per comparison)") +
  ggtitle("Top 5 DE genes per comparison (dotplot)")

print(p_dot)

# --- 8: Pathway Analysis ---

# --- 8.1 Differential Expression Analysis on groups ---

clean_DE <- function(df) {
  df <- as.data.frame(df)
  
  # store gene names
  df$gene <- rownames(df)
  rownames(df) <- NULL
  
  # convert ALL non-gene columns into numeric
  for (col in colnames(df)) {
    if (col == "gene") next
    
    # unlist if needed
    df[[col]] <- unlist(df[[col]])
    
    # convert factors to character
    if (is.factor(df[[col]])) df[[col]] <- as.character(df[[col]])
    
    # convert character to numeric when appropriate
    if (is.character(df[[col]])) {
      suppressWarnings(num <- as.numeric(df[[col]]))
      if (!all(is.na(num))) df[[col]] <- num
    }
    
    # convert Rle to numeric
    if (inherits(df[[col]], "Rle")) df[[col]] <- as.numeric(as.vector(df[[col]]))
  }
  
  # final safety cast
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  return(df)
}

DefaultAssay(integrated) <- "RNA"
Idents(integrated) <- integrated$manual_celltype

de_t_vs_b <- FindMarkers(
  integrated,
  ident.1 = "T cell (T)",
  ident.2 = "B cell (B)",
  logfc.threshold = 0,
  min.pct = 0.1,
  test.use = "wilcox"
)
de_t_vs_b <- clean_DE(de_t_vs_b)

# Extract top 5 DEGs (sorted by adjusted p-value)
top5_t_vs_b <- de_t_vs_b %>%
  arrange(p_val_adj) %>% 
  head(5) %>% 
  select(gene, avg_log2FC, p_val_adj)

print(top5_t_vs_b)

# --- 8.2 Pathway analysis on groups ---

library(dplyr)
library(enrichR)
library(ggplot2)

# 1) Select significant DE genes (T vs B)
sig_genes <- de_t_vs_b %>%
  filter(p_val_adj < 0.05) %>%
  pull(gene)

# 2) Connect to Enrichr database
enrichR::setEnrichrSite("Enrichr")

# 3) Run GO Biological Process enrichment
results <- enrichr(
  sig_genes,
  databases = c("GO_Biological_Process_2023")
)

# 4) View top enriched pathways (top 10)
head(results[["GO_Biological_Process_2023"]][, c("Term", "Adjusted.P.value", "P.value", "Overlap", "Combined.Score")], 10)

# 5) Prepare top pathways for plotting
bp <- results[["GO_Biological_Process_2023"]] %>%
  arrange(Adjusted.P.value) %>%
  head(10)

# 6) Plot top GO Biological Processes (barplot)
ggplot(bp, aes(x = reorder(Term, -Combined.Score), y = Combined.Score)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  ylab("Combined Score") +
  xlab("") +
  ggtitle("Top GO Biological Processes (T cells vs B cells)")

saveRDS(data, file ="scb.R")





























