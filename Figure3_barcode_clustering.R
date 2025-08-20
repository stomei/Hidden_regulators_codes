library(RColorBrewer)
#Set colour palettes and plot themes
brewer.pal(8, "Set1")
pop_col= c(cDC1="#E41A1C", cDC2= "#377EB8", pDC= "#4DAF4A", Mye= "#984EA3", Ery= "#FF7F00", Mast= "#FFFF33", Lymph= "#A65628", Rest= "#F781BF")

THEME1=theme(text = element_text(size = 20, colour = "black"), 
             plot.title = element_text(size = 11, face = "bold"),
             axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), 
             axis.text.y = element_text(colour = "black"), 
             axis.line = element_line(colour = "black"),
             panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             strip.background = element_blank(),
             legend.position = "right"
)

#Load dataset
all_barcodes.cpm <-read.table("ST223_all_barcodes.cpm.norm.txt", header = TRUE)
all_barcodes.cell <-read.table("ST223_all_barcodes.cell.norm.txt", header = TRUE)
common.cpm <-read.table("ST223_common.barcodes.cpm.txt", header = TRUE)
common.cell <-read.table("ST223_common.barcodes.cell.txt", header = TRUE)

#Merge plates (average values)
Samples=unique(substr(names(common.cell[,1:80]),1,nchar(names(common.cell[,1:80]))-3))
avg=data.frame(row.names=row.names(common.cell))
new.names=vector()
for (s in Samples)
{
  e=as.data.frame(common.cell[,grep(s,names(common.cell))])
  if(ncol(e)<2)next;
  new.names=c(new.names,substr(names(e)[1],1,nchar(names(e)[1])-3))
  avg=cbind(avg,rowSums(e)/2)
}

names(avg)= new.names
avg.common.cell.barcodes= avg

#run kmeans clustering
#choose numbers of centers that increases cluster compactness
library(skmeans)
clusters <- log2(avg.common.cell.barcodes[,1:40]+1)
clusters <- as.matrix(clusters)
#chose the number of centre
wssplot <- function(data, nc=15, seed=123){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of groups",
       ylab="Sum of squares within a group")}

wssplot(clusters, nc = 20) #Choose 10

clusters <- skmeans(clusters, 10)
avg.common.cell.barcodes$cluster=clusters$cluster

#Add cluster to Seurat object
ST223.rna.singlets@meta.data$barcode_donor <- paste(ST223.rna.singlets@meta.data$donor_id, ST223.rna.singlets@meta.data$barcode, sep = "_")
ST223.rna.singlets@meta.data$cell <- rownames(ST223.rna.singlets@meta.data)
all(rownames(avg.common.cell.barcodes)== rownames(common.cell))
avg.common.cell.barcodes$patient= common.cell$patient
avg.common.cell.barcodes$barcode <- substr(row.names(avg.common.cell.barcodes),4,18)
avg.common.cell.barcodes$barcode <- gsub(" ", "", avg.common.cell.barcodes$barcode)
avg.common.cell.barcodes$barcode_donor <- paste("Patient", avg.common.cell.barcodes$patient, avg.common.cell.barcodes$barcode, sep="_")
add_cluster <- select(avg.common.cell.barcodes, cluster, barcode_donor)
add_cluster <- add_cluster[!duplicated(add_cluster$barcode_donor), ]
ST223.rna.singlets@meta.data <- left_join(ST223.rna.singlets@meta.data, add_cluster, by= "barcode_donor")
rownames(ST223.rna.singlets@meta.data) <- ST223.rna.singlets@meta.data$cell

Idents(ST223.rna.singlets) <- "cluster"
DimPlot(ST223.rna.singlets, reduction = "umap")

ST223.rna.singlets.clusters <- subset(ST223.rna.singlets, idents = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
jpeg("umap_bulk_clustering.jpeg", width = 2500, height = 2000, res = 300)
DimPlot(ST223.rna.singlets.clusters, reduction = "umap")
dev.off()

ST223.rna.singlets.clusters@meta.data <- select(ST223.rna.singlets.clusters@meta.data, !seurat_clusters)
names(ST223.rna.singlets.clusters@meta.data$cluster) <- gsub("cluster", "seurat_clusters", names(ST223.rna.singlets.clusters@meta.data$cluster))
All_markers.cluster <- FindAllMarkers(ST223.rna.singlets.clusters, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

All_markers.cluster %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5_bulk

All_markers.cluster %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10_bulk

jpeg("barcode_clutsering_gene_top10_heatmap.jpeg", width = 4500, height = 3500, res = 300)
DoHeatmap(ST223.rna.singlets.clusters, features = top10_bulk$gene) + NoLegend()
dev.off()

#Find cluster based on transcriptome to compare
ST223.rna.singlets.clusters <- FindNeighbors(ST223.rna.singlets.clusters, dims = 1:20)
ST223.rna.singlets.clusters <- FindClusters(ST223.rna.singlets.clusters, resolution = 0.9)
jpeg("umap_seurat_clustering.jpeg", width = 2500, height = 2000, res = 300)
DimPlot(ST223.rna.singlets.clusters, reduction = "umap")
dev.off()

Seurat_markers.cluster <- FindAllMarkers(ST223.rna.singlets.clusters, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

Seurat_markers.cluster %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5_seurat

Seurat_markers.cluster %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10_seurat

jpeg("Seurat_clutsering_gene_heatmap.jpeg", width = 5000, height = 3500, res = 300)
DoHeatmap(ST223.rna.singlets.clusters, features = top5_seurat$gene) + NoLegend()
dev.off()

#Compare genes found with seurat clustering and bulk clustering
top10_bulk.1 <- top10_bulk[!duplicated(top10_bulk$gene), ]
top10_seurat.1 <- top10_seurat[!duplicated(top10_seurat$gene), ]
top10_bulk.1 <- top10_bulk.1[,"gene"]
top10_bulk.1$bulk <- 1
top10_seurat.1 <- top10_seurat.1[,"gene"]
top10_seurat.1$seurat <- 1
genes_comparison= full_join(top10_bulk.1, top10_seurat.1, by = "gene")

All_markers.cluster.1 <- All_markers.cluster[!duplicated(All_markers.cluster$gene), ]
Seurat_markers.cluster.1 <- Seurat_markers.cluster[!duplicated(Seurat_markers.cluster$gene), ]
All_markers.cluster.1 <- select(All_markers.cluster.1, gene)
All_markers.cluster.1$bulk <- 1
Seurat_markers.cluster.1 <- select(Seurat_markers.cluster.1, gene)
Seurat_markers.cluster.1$seurat <- 1
genes_comparison= full_join(All_markers.cluster.1, Seurat_markers.cluster.1, by = "gene")

genes_comparison[is.na(genes_comparison)]=0
genes_comparison <- as.data.frame(apply(ifelse(genes_comparison[,2:3]>0,rownames(genes_comparison),0),2,paste)) 
genes_comparison <- list('bulk'= genes_comparison$bulk,'seurat'= genes_comparison$seurat)
genes_comparison <- lapply(genes_comparison,function(x) x[x!=0])

venn.diagram(x=genes_comparison, fill= c("blue", "orange"), filename = 'cluster_comparison_venn_diagramm_all.png',
             output=TRUE)

#umap of barcodes cluster
dt <- avg.common.cpm.barcodes 
dt[,1:40] <- log2(dt[,1:40]+1)

u <- umap(dt[,1:40],metric="cosine", min_dist=0.3, n_neighbors=80, spread=15)
dt$umap1 = u$layout[,1]
dt$umap2 = u$layout[,2]
avg.common.cell.barcodes$umap1= u$layout[,1]
avg.common.cell.barcodes$umap2= u$layout[,2]

ggplot(avg.common.cell.barcodes, aes(x= umap1, y=umap2, colour = as.factor(cluster))) + 
  geom_point(size=0.5)+
  THEME1

#line plot 
line_plot <- avg.common.cell.barcodes
line_plot$barcode <- rownames(line_plot)
line_plot <-  melt(line_plot,id.vars=c("patient", "cluster", "umap1", "umap2", "barcode", "barcode_donor"))
line_plot= line_plot %>% dplyr::mutate (day = case_when(grepl("D7", variable) ~ "D7",
                                                          grepl("D10", variable)  ~ "D10",
                                                          grepl("D14", variable)  ~ "D14",
                                                          grepl("D17", variable)  ~ "D17",
                                                          grepl("D21", variable)  ~ "D21"))
line_plot$day = factor(line_plot$day, levels = c("D7","D10","D14", "D17", "D21"))
line_plot$variable <- gsub("D7_", "", line_plot$variable)
line_plot$variable <- gsub("D10_", "", line_plot$variable)
line_plot$variable <- gsub("D14_", "", line_plot$variable)
line_plot$variable <- gsub("D17_", "", line_plot$variable)
line_plot$variable <- gsub("D21_", "", line_plot$variable)
line_plot <- filter(line_plot, value >0)

line_plot <- line_plot  %>% 
  group_by(cluster, day, variable) %>%
  summarise(mean(value))
names(line_plot)[4] <- "value"

pdf(file = "line_plot_cell_common_10c.pdf", width = 6, height = 6)
cluster= unique(line_plot$cluster)
for (c in cluster)
{
  dtm = subset(line_plot, line_plot$cluster==c)
  dtm$day = factor(dtm$day, levels = c("D7","D10","D14", "D17", "D21"))
  dtm$variable = factor(dtm$variable, levels = c("cDC1","cDC2","pDC", "Mye", "Mast", "Ery", "Lymph", "Rest"))
  print(
    ggplot(dtm, aes(x=day, y=value, group = variable, colour= variable)) + 
      geom_line(linewidth = 2)+ #geom_area(alpha=0.8 , linewidth=0.1, colour="black")+
      scale_colour_manual(values = pop_col)+
      #facet_wrap(~variable)+
      THEME1 +
      ggtitle(c)
  )
}
dev.off()


FeaturePlot(ST223.rna.singlets.clusters, features = c("MPO", "LTB", "FLT3", "HBB", "SPINK2"))


write.table(avg.common.cell.barcodes, file="ST223_avg.common.cell.barcodes.txt", sep="\t", row.names=TRUE, col.names=TRUE)
