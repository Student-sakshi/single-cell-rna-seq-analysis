# ============================================================
# Human COVID-19 Immune scRNA-seq Analysis
# Single-cell transcriptomic analysis of immune-cell populations from COVID-19 and healthy samples
# Author: Sakshi Parate - M.Sc. Bioinformatics, Saarland University
# ============================================================

# ============================================================
# 1. Load Packages
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(Seurat); library(patchwork); library(DoubletFinder)
  library(SingleR); library(enrichR); library(SingleCellExperiment)
  library(celldex); library(tidyverse); library(Matrix); library(ggplot2)
  library(tidyr); library(igraph); library(leiden)
})

# ============================================================
# 2. Load Data
# ============================================================

# Input RDS files should be placed in the "data" directory.
data_dir <- "data"

sample1 <- readRDS(file.path(data_dir, "GSM4557329_GSM4557329_556_cell.counts.matrices.rds"))
sample2 <- readRDS(file.path(data_dir, "GSM4557330_GSM4557330_557_cell.counts.matrices.rds"))
sample3 <- readRDS(file.path(data_dir, "GSM4557331_GSM4557331_558_cell.counts.matrices.rds"))
sample4 <- readRDS(file.path(data_dir, "GSM4557337_GSM4557337_HIP043_cell.counts.matrices.rds"))

# ============================================================
# 3. Create Seurat Objects
# ============================================================

seurat1 <- CreateSeuratObject(counts = sample1, project = "covid_556", min.cells = 3, min.features = 200)
seurat2 <- CreateSeuratObject(counts = sample2, project = "covid_557", min.cells = 3, min.features = 200)
seurat3 <- CreateSeuratObject(counts = sample3, project = "covid_558", min.cells = 3, min.features = 200)
seurat4 <- CreateSeuratObject(counts = sample4, project = "HIP043", min.cells = 3, min.features = 200)

# ============================================================
# 4. Add Sample Metadata
# ============================================================

sample_metadata <- data.frame(
  Sample = c("covid_556", "covid_557", "covid_558", "HIP043"),
  Donor = c("C2", "C3", "C4", "H4"),
  Replicate = c("T1", "T1", "T1", "T1"),
  Sex = c("F", "M", "M", "F")
)
print(sample_metadata)

seurat1@meta.data$Donor <- "C2"; seurat1@meta.data$Replicate <- "T1"; seurat1@meta.data$Sex <- "F"
seurat2@meta.data$Donor <- "C3"; seurat2@meta.data$Replicate <- "T1"; seurat2@meta.data$Sex <- "M"
seurat3@meta.data$Donor <- "C4"; seurat3@meta.data$Replicate <- "T1"; seurat3@meta.data$Sex <- "M"
seurat4@meta.data$Donor <- "H4"; seurat4@meta.data$Replicate <- "T1"; seurat4@meta.data$Sex <- "F"

# Number of cells
ncol(seurat1); ncol(seurat2); ncol(seurat3); ncol(seurat4)

# Number of genes
nrow(seurat1); nrow(seurat2); nrow(seurat3); nrow(seurat4)

# Metadata information
colnames(seurat1@meta.data); head(seurat1@meta.data)

# ============================================================
# 5. Preprocessing
# ============================================================

# ------------------------------------------------------------
# 5.1 Quality Control
# ------------------------------------------------------------

seurat1[["percent.mt"]] <- PercentageFeatureSet(seurat1, pattern = "^MT-")
seurat2[["percent.mt"]] <- PercentageFeatureSet(seurat2, pattern = "^MT-")
seurat3[["percent.mt"]] <- PercentageFeatureSet(seurat3, pattern = "^MT-")
seurat4[["percent.mt"]] <- PercentageFeatureSet(seurat4, pattern = "^MT-")

seurat1 <- subset(seurat1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
seurat2 <- subset(seurat2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
seurat3 <- subset(seurat3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
seurat4 <- subset(seurat4, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# ------------------------------------------------------------
# 5.2 Doublet Detection
# ------------------------------------------------------------

seurat_list <- list(seurat1 = seurat1, seurat2 = seurat2, seurat3 = seurat3, seurat4 = seurat4)
df_results <- list()
expected_rate <- 0.06   # 6% expected doublet rate

for (s in names(seurat_list)) {
  message("---- Processing ", s, " ----")
  obj <- seurat_list[[s]]

  # Preprocessing required for DoubletFinder
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = 20)
  obj <- FindNeighbors(obj, dims = 1:10)
  obj <- FindClusters(obj, resolution = 0.5)

  # Parameter sweep to find optimal pK
  sweep.res <- paramSweep(obj, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res)
  pk_stats <- find.pK(sweep.stats)
  best.pK <- as.numeric(as.character(pk_stats$pK[which.max(pk_stats$BCmetric)]))
  message("Best pK for ", s, ": ", best.pK)

  # Estimate expected doublets
  nCells <- ncol(obj)
  nExp <- round(nCells * expected_rate)

  # Adjust for homotypic doublets
  homotypic.prop <- modelHomotypic(obj$seurat_clusters)
  nExp.adj <- round(nExp * (1 - homotypic.prop))

  # Run DoubletFinder
  obj <- doubletFinder(obj, PCs = 1:10, pN = 0.25, pK = best.pK, nExp = nExp.adj)

  # Identify DoubletFinder classification column
  df_col <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
  df_col <- tail(df_col, 1)
  message("DF column: ", df_col)

  # Remove doublets
  obj_clean <- obj[, obj@meta.data[[df_col]] == "Singlet"]
  message("Remaining singlet cells: ", ncol(obj_clean))

  seurat_list[[s]] <- obj_clean

  df_results[[s]] <- list(
    pK = best.pK,
    nExp = nExp,
    nExp_adjusted = nExp.adj,
    DF_column = df_col,
    cells_after = ncol(obj_clean)
  )
}

seurat1 <- seurat_list$seurat1
seurat2 <- seurat_list$seurat2
seurat3 <- seurat_list$seurat3
seurat4 <- seurat_list$seurat4
print(df_results)

# ------------------------------------------------------------
# 5.3 Normalization and Feature Selection
# ------------------------------------------------------------

seurat_list <- list(seurat1 = seurat1, seurat2 = seurat2, seurat3 = seurat3, seurat4 = seurat4)

for (s in names(seurat_list)) {
  obj <- seurat_list[[s]]
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  seurat_list[[s]] <- obj
}

seurat1 <- seurat_list$seurat1
seurat2 <- seurat_list$seurat2
seurat3 <- seurat_list$seurat3
seurat4 <- seurat_list$seurat4

# ============================================================
# 6. Batch Correction and Sample Integration
# ============================================================

# ------------------------------------------------------------
# 6.1 Merging Without Batch Correction
# ------------------------------------------------------------

merged_noBC <- merge(
  seurat1, y = list(seurat2, seurat3, seurat4),
  add.cell.ids = c("556", "557", "558", "HIP043"),
  project = "PBMC_noBatchCorrection"
)

merged_noBC <- NormalizeData(merged_noBC)
merged_noBC <- FindVariableFeatures(merged_noBC)
merged_noBC <- ScaleData(merged_noBC)
merged_noBC <- RunPCA(merged_noBC, npcs = 30)
merged_noBC <- RunUMAP(merged_noBC, dims = 1:20)

DimPlot(merged_noBC, group.by = "orig.ident", label = TRUE) +
  ggtitle("A: UMAP WITHOUT Batch Correction")

# ------------------------------------------------------------
# 6.2 Seurat Data Integration
# ------------------------------------------------------------

seurat_list <- list(seurat1, seurat2, seurat3, seurat4)

for (i in 1:length(seurat_list)) {
  seurat_list[[i]] <- NormalizeData(seurat_list[[i]])
  seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]], nfeatures = 2000)
}

features <- SelectIntegrationFeatures(seurat_list)

anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  anchor.features = features,
  dims = 1:30
)

integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, npcs = 30)
integrated <- RunUMAP(integrated, dims = 1:20)
integrated <- FindNeighbors(integrated, dims = 1:20)
integrated <- FindClusters(integrated, resolution = 0.5)

DimPlot(integrated, group.by = "orig.ident") +
  ggtitle("UMAP WITH Batch Correction (LogNormalize Integration)")

# ============================================================
# 7. Dimensionality Reduction
# ============================================================

# ------------------------------------------------------------
# 7.1 PCA and UMAP
# ------------------------------------------------------------

integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)

DimPlot(integrated, reduction = "pca", group.by = "orig.ident") +
  ggtitle("PCA Colored by Sample")

ElbowPlot(integrated, ndims = 50) +
  ggtitle("Elbow Plot for Selecting Number of PCs")

pcs_to_use <- 1:20

integrated <- RunUMAP(integrated, dims = pcs_to_use)

DimPlot(integrated, reduction = "umap", group.by = "orig.ident") +
  ggtitle("UMAP after PCA")

# ------------------------------------------------------------
# 7.2 Clustering
# ------------------------------------------------------------

integrated <- FindNeighbors(integrated, dims = pcs_to_use)

# Louvain
integrated <- FindClusters(integrated, resolution = 0.5, algorithm = 1)
Idents(integrated) <- "seurat_clusters"

DimPlot(integrated, reduction = "umap", label = TRUE) +
  ggtitle("Louvain Clustering")

# Leiden
integrated <- FindNeighbors(integrated, dims = 1:20)
names(integrated@graphs)

g <- graph_from_adjacency_matrix(
  integrated@graphs$integrated_snn,
  mode = "undirected",
  diag = FALSE,
  weighted = TRUE
)

V(g)$name <- as.character(seq_len(vcount(g)))

leiden_clusters <- leiden(g, resolution_parameter = 0.5)
integrated$leiden <- as.factor(leiden_clusters)
Idents(integrated) <- "leiden"

DimPlot(integrated, reduction = "umap", label = TRUE) +
  ggtitle("Leiden Clustering")

# ============================================================
# 8. Cell-Type Annotation
# ============================================================

# ------------------------------------------------------------
# 8.1 Automatic Cell-Type Annotation
# ------------------------------------------------------------

integrated[["integrated"]]

int_data <- GetAssayData(integrated, assay = "integrated", slot = "data")
genes <- rownames(int_data)
rna_counts <- round(exp(int_data) - 1)
rna_logdata <- int_data

sce <- SingleCellExperiment(
  assays = list(counts = rna_counts, logcounts = rna_logdata),
  colData = integrated@meta.data
)

ref <- celldex::HumanPrimaryCellAtlasData()
assayNames(ref)
ref_data <- assay(ref, "logcounts")
ref_labels <- ref$label.main
test_data <- logcounts(sce)

pred <- SingleR(
  sc_data = test_data,
  ref_data = ref_data,
  types = ref_labels,
  fine.tune = FALSE
)

integrated$SingleR_label <- pred$labels

DimPlot(integrated, reduction = "umap", group.by = "SingleR_label", label = TRUE)

# ------------------------------------------------------------
# 8.2 Manual Cell-Type Annotation
# ------------------------------------------------------------

marker.list <- list(
  HSC = c("CD34", "CD38", "SCA1", "KIT"),
  LMPP = c("CD38", "CD52", "CSF3R", "CA1", "KIT", "CD34", "FLK2"),
  CLP = c("IL7R"),
  GMP = c("ELANE"),
  CMP = c("IL3", "CSF2", "CSF1"),
  B = c("CD19", "CD20", "CD38"),
  PreB = c("CD19", "CD34"),
  Plasma = c("SDC1", "IGHA1", "IGLC1", "MZB1", "JCHAIN"),
  T = c("CD3D"),
  CD8 = c("CD3D", "CD3E", "CD8A", "CD8B"),
  CD4 = c("CD3D", "CD3E", "CD4"),
  NK = c("FCGR3A", "NCAM1", "NKG7", "KLRB1"),
  Ery = c("GATA1", "HBB", "HBA1", "HBA2"),
  pDC = c("IRF8", "IRF4", "IRF7"),
  cDC = c("CD1C", "CD207", "ITGAM", "NOTCH2", "SIRPA"),
  Mono14 = c("CD14", "CCL3", "CCL4", "IL1B"),
  Mono16 = c("FCGR3A", "CD68", "S100A12"),
  Basophils = c("GATA2")
)

marker.list_lower <- lapply(marker.list, function(x) tolower(x))

expr_mat <- GetAssayData(integrated, assay = "RNA", slot = "data")
gene_names <- rownames(expr_mat)
gene_names_lower <- tolower(gene_names)

marker_hits_per_type <- sapply(
  marker.list_lower,
  function(m) sum(m %in% gene_names_lower)
)

print("marker hits per marker-set:")
print(marker_hits_per_type)

Idents(integrated) <- integrated$seurat_clusters
cluster_ids <- levels(Idents(integrated))
cluster_to_celltype <- setNames(rep(NA_character_, length(cluster_ids)), cluster_ids)

for (cl in cluster_ids) {
  cells <- WhichCells(integrated, idents = cl)
  if (length(cells) == 0) next

  avg <- Matrix::rowMeans(expr_mat[, cells, drop = FALSE])
  names(avg) <- gene_names

  scores <- sapply(marker.list_lower, function(marker_genes_lower) {
    hits_idx <- which(gene_names_lower %in% marker_genes_lower)
    if (length(hits_idx) == 0) return(0)
    mean(avg[hits_idx], na.rm = TRUE)
  })

  best <- names(scores)[which.max(scores)]
  cluster_to_celltype[cl] <- best
}

print("cluster -> chosen marker set (abbrev):")
print(cluster_to_celltype)

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

cluster_pretty <- sapply(
  cluster_to_celltype,
  function(abbrev) {
    if (is.na(abbrev)) return("Unknown")
    lab <- celltype.pretty[[abbrev]]
    if (is.null(lab)) {
      return(paste0(toupper(substring(abbrev, 1, 1)), substring(abbrev, 2), " (", abbrev, ")"))
    } else {
      return(lab)
    }
  },
  USE.NAMES = FALSE
)

manual_labels_map <- setNames(cluster_pretty, cluster_ids)

integrated$manual_celltype <- unname(
  manual_labels_map[as.character(Idents(integrated))]
)

message("Manual annotation counts:")
print(table(integrated$manual_celltype, useNA = "ifany"))

DimPlot(
  integrated,
  reduction = "umap",
  group.by = "manual_celltype",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Manual Cell Type Annotation (Table 2 markers)")

if ("SingleR_label" %in% colnames(integrated@meta.data)) {
  p1 <- DimPlot(integrated, group.by = "SingleR_label", label = TRUE) + ggtitle("Automatic")
  p2 <- DimPlot(integrated, group.by = "manual_celltype", label = TRUE) + ggtitle("Manual")
  print(p1 + p2)
}

markers_to_plot <- c("CD3D", "MS4A1", "NKG7")

plots <- VlnPlot(
  integrated,
  features = markers_to_plot,
  group.by = "seurat_clusters",
  pt.size = 0,
  combine = FALSE
)

wrap_plots(plots, ncol = 3) +
  plot_annotation(title = "Expression of 3 marker genes")

FeaturePlot(integrated, features = markers_to_plot, reduction = "umap")

Idents(integrated) <- integrated$manual_celltype
integrated$celltype_merged <- Idents(integrated)

message("Finished manual annotation and merging.")

# ------------------------------------------------------------
# 8.3 Cell-Type Proportions
# ------------------------------------------------------------

Idents(integrated) <- integrated$manual_celltype
integrated$sample <- integrated$orig.ident

prop_df <- integrated@meta.data %>%
  group_by(sample, manual_celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count))

print(prop_df)

ggplot(prop_df, aes(x = sample, y = proportion, fill = manual_celltype)) +
  geom_bar(stat = "identity") +
  ylab("Proportion") +
  xlab("Sample") +
  ggtitle("Cell-Type Proportions per Sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================
# 9. Differential Expression Analysis
# ============================================================

# ------------------------------------------------------------
# 9.1 Differential Expression on Cell Types
# ------------------------------------------------------------

raw_counts <- GetAssayData(integrated, assay = "RNA", slot = "counts")

dups <- duplicated(rownames(raw_counts))
message("Duplicate genes removed: ", sum(dups))

raw_counts_clean <- raw_counts[!dups, ]
newRNA <- CreateAssayObject(counts = raw_counts_clean)
integrated[["RNA"]] <- newRNA

DefaultAssay(integrated) <- "RNA"

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

if (!dir.exists("results")) dir.create("results")

write.csv(de_b_vs_t, file.path("results", "DE_B_vs_T.csv"), row.names = FALSE)
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

write.csv(de_t_vs_mono, file.path("results", "DE_T_vs_Mono16.csv"), row.names = FALSE)
message("T vs Mono16 DE complete.")

# ------------------------------------------------------------
# 9.2 Volcano Plots
# ------------------------------------------------------------

volcano_b_t <- ggplot(de_b_vs_t, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-0.5, 0.5), color = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed") +
  ggtitle("Volcano Plot: B cells vs T cells") +
  xlab("log2 Fold Change") +
  ylab("-log10 adjusted p-value") +
  theme_bw()

print(volcano_b_t)

volcano_t_mono <- ggplot(de_t_vs_mono, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-0.5, 0.5), color = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed") +
  ggtitle("Volcano Plot: T cells vs Monocyte (CD16)") +
  xlab("log2 Fold Change") +
  ylab("-log10 adjusted p-value") +
  theme_bw()

print(volcano_t_mono)

# ------------------------------------------------------------
# 9.3 COVID-19 vs Healthy Differential Expression
# ------------------------------------------------------------

integrated$condition <- ifelse(integrated$orig.ident == "HIP043", "Healthy", "COVID")
Idents(integrated) <- integrated$condition
levels(Idents(integrated))

de_covid_vs_healthy <- FindMarkers(
  integrated,
  ident.1 = "COVID",
  ident.2 = "Healthy",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"
)

de_covid_vs_healthy$gene <- rownames(de_covid_vs_healthy)

write.csv(
  de_covid_vs_healthy,
  file.path("results", "DE_COVID_vs_Healthy.csv"),
  row.names = FALSE
)

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

print(top_up)
print(top_down)

print(head(de_covid_vs_healthy[order(de_covid_vs_healthy$p_val_adj), ], 10))

# ------------------------------------------------------------
# 9.4 Top Differentially Expressed Genes
# ------------------------------------------------------------

clean_DE <- function(df) {
  df <- as.data.frame(df)
  df$gene <- rownames(df)
  rownames(df) <- NULL

  for (col in colnames(df)) {
    if (inherits(df[[col]], "Rle")) df[[col]] <- as.vector(df[[col]])
    if (is.list(df[[col]])) df[[col]] <- unlist(df[[col]])
    if (is.factor(df[[col]])) df[[col]] <- as.character(df[[col]])

    if (is.character(df[[col]])) {
      suppressWarnings(num <- as.numeric(df[[col]]))
      if (!all(is.na(num))) df[[col]] <- num
    }
  }

  df <- as.data.frame(df, stringsAsFactors = FALSE)
  return(df)
}

de_b_vs_t <- clean_DE(de_b_vs_t)
de_t_vs_mono <- clean_DE(de_t_vs_mono)

if (!"p_val_adj" %in% colnames(de_b_vs_t) & "p_val" %in% colnames(de_b_vs_t))
  de_b_vs_t$p_val_adj <- de_b_vs_t$p_val

if (!"p_val_adj" %in% colnames(de_t_vs_mono) & "p_val" %in% colnames(de_t_vs_mono))
  de_t_vs_mono$p_val_adj <- de_t_vs_mono$p_val

de_b_vs_t <- de_b_vs_t %>%
  mutate(
    avg_log2FC = as.numeric(avg_log2FC),
    p_val_adj = as.numeric(p_val_adj),
    neglog10_padj = -log10(pmax(p_val_adj, 1e-300)),
    enriched_in = ifelse(avg_log2FC > 0, "B cell (B)", "T cell (T)")
  )

de_t_vs_mono <- de_t_vs_mono %>%
  mutate(
    avg_log2FC = as.numeric(avg_log2FC),
    p_val_adj = as.numeric(p_val_adj),
    neglog10_padj = -log10(pmax(p_val_adj, 1e-300)),
    enriched_in = ifelse(avg_log2FC > 0, "T cell (T)", "Monocyte CD16 (Mono16)")
  )

top5_b_t <- de_b_vs_t %>%
  arrange(p_val_adj) %>%
  slice_head(n = 5) %>%
  mutate(comparison = "B_vs_T")

top5_t_mono <- de_t_vs_mono %>%
  arrange(p_val_adj) %>%
  slice_head(n = 5) %>%
  mutate(comparison = "T_vs_Mono16")

top5_all <- bind_rows(top5_b_t, top5_t_mono)

top5_all <- top5_all %>%
  group_by(comparison) %>%
  mutate(gene = factor(gene, levels = rev(gene))) %>%
  ungroup()

p_dot <- ggplot(
  top5_all,
  aes(x = enriched_in, y = gene, size = neglog10_padj, color = avg_log2FC)
) +
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

# ============================================================
# 10. Pathway Analysis
# ============================================================

# ------------------------------------------------------------
# 10.1 Differential Expression: T Cells vs B Cells
# ------------------------------------------------------------

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

top5_t_vs_b <- de_t_vs_b %>%
  arrange(p_val_adj) %>%
  head(5) %>%
  select(gene, avg_log2FC, p_val_adj)

print(top5_t_vs_b)

# ------------------------------------------------------------
# 10.2 GO Biological Process Enrichment
# ------------------------------------------------------------

enrichR::setEnrichrSite("Enrichr")

sig_genes <- de_t_vs_b %>%
  filter(p_val_adj < 0.05) %>%
  pull(gene)

results <- enrichr(
  sig_genes,
  databases = c("GO_Biological_Process_2023")
)

go_results <- results[["GO_Biological_Process_2023"]] %>%
  select(Term, Adjusted.P.value, P.value, Overlap, Combined.Score)

print(head(go_results, 10))

write.csv(
  go_results,
  file.path("results", "GO_Biological_Process_Enrichment.csv"),
  row.names = FALSE
)

bp <- results[["GO_Biological_Process_2023"]] %>%
  arrange(Adjusted.P.value) %>%
  head(10)

ggplot(
  bp,
  aes(x = reorder(Term, -Combined.Score), y = Combined.Score)
) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  ylab("Combined Score") +
  xlab("") +
  ggtitle("Top GO Biological Processes (T cells vs B cells)")

# ============================================================
# End of Analysis
# ============================================================
