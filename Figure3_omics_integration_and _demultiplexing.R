#install packages
install.packages("spgs")

#load libraries
library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)
library(SingleR)
library(ComplexHeatmap)
library(DoMultiBarHeatmap)
library(rlang)
library(readr)
library(maditr)
library(reshape2)
library(spgs)
library(dplyr)
library(tidyr)

#RNA alanysis only
#Load dataset from mgi
ST223.mgi <- Read10X(data.dir = "/wehisan/general/user_managed/grpu_naik.s_2/2023_Sequencing_Runs/ST223/Analsyis/ST223_GEX/outs/filtered_feature_bc_matrix/")
ST223.rna <- CreateSeuratObject(counts = ST223.mgi, project = "ST223")

#Add SPLINTR to RNA dataset
SPLINTR.1 <-read.table("/stornext/General/data/academic/lab_naik/Sara_Tomei/R_analysis/10X_Analysis/ST223/V350096710_S1_L001_splintr_counts.csv", header = FALSE, sep = ",")
SPLINTR.2 <-read.table("/stornext/General/data/academic/lab_naik/Sara_Tomei/R_analysis/10X_Analysis/ST223/V350096710_S1_L002_splintr_counts.csv", header = FALSE, sep = ",")
meta.data <-read.table("/stornext/General/data/academic/lab_naik/Sara_Tomei/R_analysis/10X_Analysis/ST223/meta.data.txt", header = TRUE, sep = "\t")
SPLINTR <- rbind(SPLINTR.1, SPLINTR.2)
SPLINTR <- select(SPLINTR, V1, V2, V4)
SPLINTR <- SPLINTR[!duplicated(SPLINTR), ]
SPLINTR$V4 <- gsub("GCCCTGAT", "", SPLINTR$V4)
SPLINTR$V4 <- gsub("TATGCAAG", "", SPLINTR$V4)
SPLINTR=data.frame(SPLINTR, barcode=substr(SPLINTR$V4,1,15))
barcodes <- SPLINTR(select, V1, barcode)
barcodes <- barcodes[!duplicated(barcodes$V1), ]
barcodes$cell <- paste(barcodes$V1, "-1", sep="")

meta.data <- select(meta.data, cell, donor_id, cons_bc)

ST223.rna@meta.data$cell=rownames(ST223.rna@meta.data)
ST223.rna@meta.data <- left_join(ST223.rna@meta.data, meta.data, by="cell")

#Add demultiplex with SNP 
donors <- read.table("/wehisan/general/user_managed/grpu_naik.s_2/2023_Sequencing_Runs/ST223/Analsyis/GEX_demux/donor_ids.tsv", header = TRUE, sep ="\t")
donors <- select (donors, cell, donor_id)
ST223.rna@meta.data <- left_join(ST223.rna@meta.data, donors, by="cell")
rownames(ST223.rna@meta.data) <- ST223.rna@meta.data$cell
#ST223.rna@meta.data <- select(ST223.rna@meta.data, !cell_ID)

#Check barcodes in common with dataset
Total_barcodes <- full_join(ST218.1, ST208.1, by = "barcode")
Total_barcodes <- separate(Total_barcodes, barcode, c("barcode", "protein"))
Total_barcodes <- select(Total_barcodes, barcode, protein)
write.table(Total_barcodes, file= "Total_barcodes.txt", sep="\t", row.names=FALSE, col.names=TRUE)
Total_barcodes=data.frame(Total_barcodes, barcode=substr(Total_barcodes$barcode,1,20))
Total_barcodes=select(Total_barcodes, !barcode)
names(Total_barcodes)[2] <- "barcode"

try <- inner_join(Total_barcodes, SPLINTR, by ="barcode")

#Add ADT and HTO to RNA
mbm.ST223.adt <- Read10X(data.dir = "/wehisan/general/user_managed/grpu_naik.s_2/2023_Sequencing_Runs/ST223/Analsyis/ST223_ADT_HTO/outs/raw_feature_bc_matrix/")
ST223.adt.illumina <- Read10X(data.dir = "/wehisan/general/user_managed/grpu_naik.s_2/2023_Sequencing_Runs/ST223/Analsyis/ST223_ADT_HTO_re_run/outs/raw_feature_bc_matrix/")

HTO <- RowMergeSparseMatrices(mbm.ST223.adt, ST223.adt.illumina )
ADT <- readRDS("/stornext/General/data/academic/lab_naik/Sara_Tomei/R_analysis/10X_Analysis/ST223/ADT.rds")

hto.counts <- HTO[grepl("HTO", row.names(HTO)), ]
adt.counts <- ADT@counts

joint.bcs <- intersect(colnames(ST223.rna), colnames(adt.counts))
ST223.rna <- ST223.rna[, joint.bcs]
#hto.counts <- hto.counts[, joint.bcs]
adt.counts <- adt.counts[, joint.bcs]

adt_assay <- CreateAssayObject(counts = adt.counts)
hto_assay <- CreateAssayObject(counts = hto.counts)


# add this assay to the previously created Seurat object
ST223.rna[["ADT"]] <- adt_assay
ST223.rna[["HTO"]] <- hto_assay

# Validate that the object now contains multiple assays
Assays(ST223.rna)

#Confirm that HTo have the right name
rownames(ST223.rna@assays$HTO)

#Filter out low quality data
VlnPlot(ST223.rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ST223.rna <- subset(ST223.rna, subset = nFeature_RNA > 200 & nFeature_RNA < 7500)

#Demultiplex the sample
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
try <- subset(ST223.rna, subset = HTO-13 > 0)
ST223.rna <- NormalizeData(ST223.rna, assay = "HTO", normalization.method = "CLR")
#Use HTODemux
ST223.rna<- HTODemux(ST223.rna, assay = "HTO", positive.quantile = 0.99)
#Visualize demultiplexing results
# Global classification results
table(ST223.rna$HTO_classification.global)
rowSums(ST223.rna@assays$HTO)

#Reads count for singlets and doublets
Idents(ST223.rna) <- "HTO_classification.global"
VlnPlot(ST223.rna, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

#Visulaize HTO expression level
# Group cells based on the max HTO signal
Idents(ST223.rna) <- "HTO_maxID"
RidgePlot(ST223.rna, assay = "HTO", features = rownames(ST223.rna[["HTO"]]), ncol = 2)

# HTO heatmap
HTOHeatmap(ST223.rna, assay = "HTO")

# Extract the singlets
Idents(ST223.rna) <- "donor_id"
ST223.rna.singlets <- subset(ST223.rna, idents = "doublet", invert = TRUE)

saveRDS(ST223.rna.adt.singlets, file = "/stornext/General/data/academic/lab_naik/Sara_Tomei/R_analysis/10X_Analysis/ST223/ST223.annotated.rds")
