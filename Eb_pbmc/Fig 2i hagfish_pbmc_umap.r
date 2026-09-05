library(Seurat)
library(SCpubr)
pbmc <- readRDS("Eb_pbmc.rds")
pbmc$group_clusters <- paste0(pbmc$custom_clusters, "_", pbmc$seurat_clusters)
#color set
colours <- c(
"#FF9A76",
"#8963BA",
"#34688D",
"#6CB9A6",
"#C1719C"
)
color_plot <- colours
names(color_plot)<-unique(Idents(pbmc))
do_DimPlot(pbmc,pt.size=0.1,reduction ="umap",colors.use = color_plot) 

cell_counts <- table(Idents(pbmc))

cell_prop <- as.data.frame(cell_counts)
colnames(cell_prop) <- c("Cluster", "Cell_Number")
cell_prop$Proportion <- cell_prop$Cell_Number / sum(cell_prop$Cell_Number)

cell_prop
#       Cluster Cell_Number Proportion
#1     VLRA+/C+       20901 0.68294994
#2  Thrombocyte        6252 0.20428702
#3        VLRB+        2629 0.08590380
#4     Monocyte         389 0.01271076
#5 Erythorocyte         433 0.01414848