# install sctransform
install.packages("sctransform")
BiocManager::install("glmGamPoi")
BiocManager::install("ComplexHeatmap")
BiocManager::install("DoMultiBarHeatmap")

install.packages("class", "KernSmooth", "MASS", "Matrix", "nnet")
#load packages
library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(DoMultiBarHeatmap)

#Import ST223.rna.singlets from ST223_omics_integration_and_demultiplexing
ST223.rna.singlets <- readRDS("/stornext/General/data/academic/lab_naik/Sara_Tomei/R_analysis/10X_Analysis/ST223/ST223.rna.singlets.rds")

# Filter out poor quality cells
ST223.rna.singlets[["percent.mt"]] <- PercentageFeatureSet(ST223.rna.singlets, pattern = "^MT-")
VlnPlot(ST223.rna.singlets, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ST223.rna.singlets <- subset(ST223.rna.singlets, idents = "unassigned", invert = TRUE)
ST223.rna.singlets <- subset(ST223.rna.singlets, subset = percent.mt < 20)

#Normalise data
ST223.rna.singlets <- NormalizeData(ST223.rna.singlets) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ST223.rna.singlets <- FindNeighbors(ST223.rna.singlets, dims = 1:10)
ST223.rna.singlets <- FindClusters(ST223.rna.singlets, resolution = 0.5)
ST223.rna.singlets <- RunUMAP(ST223.rna.singlets, reduction = "pca", dims = 1:30, verbose = FALSE)

Idents(ST223.rna.singlets) <- "donor_id"
jpeg("SNP_donors.jpeg", width = 2000, height = 1500, res = 300)
DimPlot(ST223.rna.singlets, reduction = "umap")
dev.off()
FeaturePlot(ST223.rna.singlets, features = donor_id)

#Get rid of Esmaeel samples and re-do normalisation
ST223.rna.singlets <- subset(ST223.rna.singlets, idents = "donor3", invert = TRUE)
ST223.rna.singlets <- NormalizeData(ST223.rna.singlets) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ElbowPlot(ST223.rna.singlets)

ST223.rna.singlets <- FindNeighbors(ST223.rna.singlets, dims = 1:20)
ST223.rna.singlets <- FindClusters(ST223.rna.singlets, resolution = 0.2)
ST223.rna.singlets <- RunUMAP(ST223.rna.singlets, reduction = "pca", dims = 1:20, verbose = FALSE)

jpeg("umap_cluster).jpeg", width = 2000, height = 1500, res = 300)
DimPlot(ST223.rna.singlets, reduction = "umap")
dev.off()

#ST223.rna.singlets  <- SCTransform(ST223.rna.singlets , vst.flavor = "v2", verbose = FALSE) %>%
  #RunPCA(npcs = 20, verbose = FALSE) %>%
  #RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
  #FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
  #FindClusters(resolution = 0.7, verbose = FALSE)

#Find marker genes
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
All_markers <- FindAllMarkers(ST223.rna.singlets, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
All_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

jpeg("gene_heatmap.jpeg", width = 3000, height = 3000, res = 300)
DoHeatmap(ST223.rna.singlets, features = top10$gene) + NoLegend()
dev.off()

#Find which Donor is which patient based on barcodes
barcodes.1 <- left_join(donors, barcodes, by = "cell")
barcodes.1 <- barcodes.1 %>% filter(!is.na(barcode))
Tot_barcodes=read.table("/stornext/General/data/academic/lab_naik/Sara_Tomei/R_analysis/Barcode_Analysis/ST218/ST218_Total_barcode_norm_cell_for_scRNAseq.txt",header=TRUE)
Tot_barcodes=data.frame(Tot_barcodes, barcode=substr(rownames(Tot_barcodes),4,18))
Tot_barcodes$barcode=gsub(" ", "", Tot_barcodes$barcode)
Tot_barcodes=Tot_barcodes[rowSums(Tot_barcodes[,1:120])>0,]
Tot_barcodes_P=filter(Tot_barcodes, patient =="P")
Tot_barcodes_1=filter(Tot_barcodes, patient =="1")
Tot_barcodes <- Tot_barcodes[!duplicated(Tot_barcodes$barcode), ]
Donor0 <- filter(barcodes.1, donor_id== "donor0")
Donor0 <- left_join(Donor0, Tot_barcodes_1, by="barcode")
Donor0 <- Donor0 %>% filter(!is.na(patient)) #P
Donor1 <- filter(barcodes.1, donor_id== "donor1")
Donor1 <- left_join(Donor1, Tot_barcodes, by="barcode")
Donor1 <- Donor1 %>% filter(!is.na(patient)) #V
Donor2 <- filter(barcodes.1, donor_id== "donor2")
Donor2 <- left_join(Donor2, Tot_barcodes, by="barcode")
Donor2 <- Donor2 %>% filter(!is.na(patient)) #U
Donor4 <- filter(barcodes.1, donor_id== "donor4") #1
Donor4 <- left_join(Donor4, Tot_barcodes_1, by="barcode")
Donor4 <- Donor4 %>% filter(!is.na(patient))

ST223.rna.singlets@meta.data$donor_id <- gsub("donor0", "P", ST223.rna.singlets@meta.data$donor_id)
ST223.rna.singlets@meta.data$donor_id <- gsub("donor1", "V", ST223.rna.singlets@meta.data$donor_id)
ST223.rna.singlets@meta.data$donor_id <- gsub("donor2", "U", ST223.rna.singlets@meta.data$donor_id)
ST223.rna.singlets@meta.data$donor_id <- gsub("donor4", "1", ST223.rna.singlets@meta.data$donor_id)

Idents(ST223.rna.singlets) <- "donor_id"
jpeg("SNP_filtered_donors.jpeg", width = 2000, height = 1500, res = 300)
DimPlot(ST223.rna.singlets, reduction = "umap")
dev.off()

#Count barcodes shared per plate
Donor0$patient <- "P"
Donor1$patient <- "V"
Donor2$patient <- "U"
Donor4$patient <- "p1"

Patinets_barcodes <- rbind(Donor0, Donor1, Donor2, Donor4)
patients= unique(Patinets_barcodes$patient)
plate= c("PA", "PB", "PC")

d=data.frame()
for (c in patients)
{
  dtm= subset(Patinets_barcodes, Patinets_barcodes$patient==c)
  tot_barcode_num= nrow(dtm)
  b=data.frame()
  for (p in plate)
  {
    e=as.data.frame(dtm[,grep(p,names(dtm),value =TRUE )])
    e[is.na(e)]=0
    e=e[rowSums(e)>0,]
    barcode_num= nrow(e)
    percentage= nrow(e)/nrow(dtm)*100
    a=data.frame(patient= c, plate= p, barcode_num= barcode_num, tot_barcode_num= tot_barcode_num, percentage= percentage)
    b= rbind(b, a)
  }
  d = rbind(d,b)
}

Idents(ST223.rna.singlets) <- "seurat_clusters"
#ComplexHeatmap
counts.mgi <- GetAssayData(ST223.rna.singlets, assay="RNA", slot="data")
counts.mgi <- as.matrix(counts.mgi[rownames(counts.mgi) %in% top10$gene, ])
metadata.mgi <- ST223.rna.singlets[["donor_id"]]
metadata.mgi$cell <- rownames(metadata.mgi)
metadata.mgi <- metadata.mgi[rownames(metadata.mgi) %in% colnames(counts.mgi), ]

ann <- data.frame(metadata.mgi$donor_id)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col')
Heatmap(counts.mgi, name = "mat", top_annotation=colAnn)
mouse <- HeatmapAnnotation(as.factor(metadata.mgi))

anno_row = metadata.mgi$HTO_classification
pheatmap(data.matrix(asinh(counts.mgi)),
         cluster_rows = T, cluster_cols = T,annotation_col = as.factor(metadata.mgi$HTO_classification),
         color=colorRampPalette(c("grey80", "red"))(10),
         border_color = NA,
         legend = F, annotation_legend = T, show_rownames = F, fontsize = 7,
         clustering_distance_cols = "correlation",)

DoMultiBarHeatmap(object = m.bm.singlets.mgi, features = top5$gene, group.by="cluster", size = 2 ) + ggtitle("data2_mygenelist")

try <- as.data.frame(t(ST223.rna.singlets@assays$RNA@counts), patient= ST223.rna.singlets@meta.data$donor_id, clusters= ST223.rna.singlets@meta.data$seurat_clusters)
