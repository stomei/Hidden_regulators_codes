install.packages("SeuratData")
devtools::install_github('satijalab/seurat-data')
BiocManager::install("ComplexHeatmap")

#load packages
library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(DoMultiBarHeatmap)

#Import ST223.rna.adt.singlets from ST223_omics_integration_and_demultiplexing
ST223.rna.adt.singlets <- readRDS("/stornext/General/data/academic/lab_naik/Sara_Tomei/R_analysis/10X_Analysis/ST223/ST223.rna.adt.singlets.rds")

# Filter out poor quality cells
ST223.rna.adt.singlets[["percent.mt"]] <- PercentageFeatureSet(ST223.rna.adt.singlets, pattern = "^MT-")
VlnPlot(ST223.rna.adt.singlets, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ST223.rna.adt.singlets <- subset(ST223.rna.adt.singlets, subset = percent.mt < 20)

#Normalise data
DefaultAssay(ST223.rna.adt.singlets) <- 'RNA'
ST223.rna.adt.singlets <- NormalizeData(ST223.rna.adt.singlets) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(ST223.rna.adt.singlets) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the
VariableFeatures(ST223.rna.adt.singlets) <- rownames(ST223.rna.adt.singlets[["ADT"]])
ST223.rna.adt.singlets <- NormalizeData(ST223.rna.adt.singlets, normalization.method = 'CLR', margin = 2) %>%
  ScaleData() %>% RunPCA(reduction.name = 'apca')

#Find Neighbours using both RNA and ADT
ST223.rna.adt.singlets <- FindMultiModalNeighbors(
  ST223.rna.adt.singlets, reduction.list = list("pca", "apca"),
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)

VlnPlot(ST223.rna.adt.singlets, features = "RNA.weight", group.by = 'seurat_clusters', sort = TRUE, pt.size = 1) +
  NoLegend()

#ST223.rna.adt.singlets <- RunUMAP(ST223.rna.adt.singlets, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
#ST223.rna.adt.singlets <- FindClusters(ST223.rna.adt.singlets, graph.name = "wsnn", algorithm = 3, resolution = 1.1, verbose = FALSE)

DimPlot(ST223.annotated, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 4) + NoLegend()

ST223.rna.adt.singlets <- RunUMAP(ST223.rna.adt.singlets, reduction = 'pca', dims = 1:30, assay = 'RNA',
              reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
ST223.rna.adt.singlets <- RunUMAP(ST223.rna.adt.singlets, reduction = 'apca', dims = 1:18, assay = 'ADT',
              reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

DimPlot(ST223.rna.adt.singlets, reduction = 'rna.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()

DefaultAssay(ST223.rna.adt.singlets) <- 'RNA'
ST223_cluster_marker <- FindAllMarkers(ST223.rna.adt.singlets, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

ST223_markers_10 <- ST223_cluster_marker %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

jpeg("heatmap_clustering.jpeg", width = 2500, height = 6000, res = 300)
DoHeatmap(ST223.rna.adt.singlets, features = ST223_markers_10$gene) + NoLegend()
dev.off()

DefaultAssay(ST223.rna.adt.singlets) <- 'ADT'
ST223_adt_marker <- FindAllMarkers(ST223.rna.adt.singlets, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

ST223_adt_5 <- ST223_adt_marker %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

jpeg("heatmap_adt_clustering.jpeg", width = 3000, height = 5000, res = 300)
DoHeatmap(ST223.rna.adt.singlets, features = ST223_adt_5$gene) + NoLegend()
dev.off()

FeaturePlot(ST223.rna.adt.singlets, reduction = 'wnn.umap', features = c("Hu.CD22", "Hu.CD116", "Hu.CX3CR1", "Hu.CD123", "Hu.CD235"))


=# you can plot raw counts as well
VlnPlot(ST223.rna.adt.singlets, features = c( "PLD4", "PLAC8", "CCDC50", "IL3RA", 'IRF8', 'KCNE5', 'SCT'), slot = "counts", log = TRUE)
VlnPlot(ST223.rna.adt.singlets, features = c( "PRTN3", "AZU1", "PLD4", "S100A9", "PLAC8", "STAG3", "MNDA" ,"ELANE", "S100A12", "LYZ", "RETN", "MPO", "RNASE2") , slot = "counts", log = TRUE)
VlnPlot(ST223.rna.adt.singlets, features = c("ZCCHC7", 'ADA', 'IRF1', 'LTB', 'IGLL1', 'UHRF1', 'PSME1', 'CDK6', 'SEPTIN6', 'LAT2', 'ETS2', 'LSP1', 'SOX4'), slot = "counts", log = TRUE)
VlnPlot(ST223.rna.adt.singlets, features = c('IRF1', 'ZCCHC7', 'ADA', 'SEPTIN6', 'CDK6', 'IGLL1'), slot = "counts", log = TRUE)
VlnPlot(ST223.rna.adt.singlets, features = c('CA1', 'HBB', 'HBD', 'AHSP', 'KLF1', 'CAVIN2', 'ANK1', 'TUBA1B', 'HIST1H4C', 'ITGA2B', 'ICAM4', "DENND10"),slot = "counts", log = TRUE)
VlnPlot(ST223.rna.adt.singlets, features = c('LTB',"LSP1", "MPO","XIST", "ZCCHC7"),slot = "counts", log = TRUE)
VlnPlot(ST223.rna.adt.singlets, features = c('ITGA2B', "TASL", "RPS6KA4", "PLD4", 'IRF1', "SEPTIN6"),slot = "counts", log = TRUE)
VlnPlot(ST223.rna.adt.singlets, features = c('LMO4', 'IL5RA', 'CLC', 'PRG2', 'GATA2', 'HPGDS', 'HDC', 'CPA3', 'RHEX', 'CAMLG'),slot = "counts", log = TRUE)

new.cluster.ids <- c("Ery_prog", "Mk", "Ery_prog", "early_lymphoid_progenitor", "HSC", "Mk-ery_prog",
                     "basophil-mast_progenitor", "pre-cDC", "Mk", "DC", "MDP", "CMP", "CD14+Monocyte", "MPP", "MPP", "Mk", "16")
new.cluster.ids <- c("primitive-RBC", "Mk-ery_prog", "Mk", "early_lymphoid_progenitor", "HSC", "basophil-mast_progenitor", "pre-cDC", "Ery", "Mk", "CDP", "CD14+Monocyte", "CMP-GMP",  "Ery", "Ery", "Mk", "15")
names(new.cluster.ids) <- levels(ST223.rna.adt.singlets)
ST223.rna.adt.singlets <- RenameIdents(ST223.rna.adt.singlets, new.cluster.ids)
DimPlot(ST223.annotated_fate, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#reference dataset
#load reference data
InstallData("bmcite")
bm <- LoadData(ds = "bmcite")

VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>%
  ScaleData() %>% RunPCA(reduction.name = 'apca')

bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap",
              reduction.key = "wnnUMAP_", return.model = TRUE)
DimPlot(bm, group.by = "celltype.l2", reduction = "wnn.umap")

Idents(bm) <- "celltype.l2"
bm.markers <- FindAllMarkers(bm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

bm_markers_10 <- bm.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

DoHeatmap(bm, features = bm_markers_10$gene) + NoLegend()

DefaultAssay(bm) <- 'ADT'
bm.adt.markers <- FindAllMarkers(bm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

bm_adt_markers_5 <- bm.adt.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

DoHeatmap(bm, features = bm_adt_markers_5$gene) + NoLegend()

bm <- ScaleData(bm, assay = 'RNA')
bm <- RunSPCA(bm, assay = 'RNA', graph = 'wsnn')

ST223.rna.adt.singlets <- RunSPCA(ST223.rna.adt.singlets, assay = 'RNA', graph = 'wsnn')

bm <- FindNeighbors(
  object = bm,
  reduction = "spca",
  dims = 1:50,
  graph.name = "spca.annoy.neighbors",
  k.param = 50,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

#Find anchors
DefaultAssay(ST223.rna.adt.singlets) <- 'RNA'
DefaultAssay(bm) <- 'RNA'

anchors <- FindTransferAnchors(
           reference = bm,
           query = ST223.rna.adt.singlets,
           k.filter = NA,
           reference.reduction = "spca",
           reference.neighbors = "spca.annoy.neighbors",
           dims = 1:50)

ST223.rna.adt.singlets <- MapQuery(
    anchorset = anchors,
    query = ST223.rna.adt.singlets,
    reference = bm,
    refdata = list(
      celltype = "celltype.l2",
      predicted_ADT = "ADT"),
    reference.reduction = "spca",
    reduction.model = "wnn.umap")

p1 <- DimPlot(ST223.rna.adt.singlets, reduction = 'wnn.umap', group.by = 'predicted.celltype', label.size = 3)
p2 <- DimPlot(ST223.rna.adt.singlets, reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3)
p1 + p2 + plot_layout(guides = "collect")



p1

