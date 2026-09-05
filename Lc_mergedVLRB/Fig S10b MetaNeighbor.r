library(MetaNeighbor)
library(Seurat)
library(ComplexHeatmap)
library(SingleCellExperiment)
VLRB
VLRB_sub <- subset(VLRB, subset = group != "Gi")
VLRB.sce <- as.SingleCellExperiment(VLRB_sub, assay = "SCT")
global_hvgs1<-variableGenes(dat=VLRB.sce,exp_labels=VLRB.sce$group)
length(global_hvgs1)
Aurocs_matrix=MetaNeighborUS(var_genes=global_hvgs1,
dat=VLRB.sce,
study_id=VLRB.sce$group,
cell_type=VLRB.sce$seurat_clusters,
fast_version=T)
library(pheatmap)

d <- 1 - Aurocs_matrix   # AUROC → distance

pheatmap(Aurocs_matrix,
         clustering_distance_rows = as.dist(d),
         clustering_distance_cols = as.dist(d),
         clustering_method = "ward.D2",
         color = colorRampPalette(c("grey", "white", "#3288BD"))(50), #c("#3288BD", "white", "#D53E4F")
         rect_gp = gpar(col = "black"),
fontsize_row = 8,
fontsize_col = 8)#7.42*6.52 #portrait 