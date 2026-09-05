library(SCpubr)
library(Seurat)

pbmc <- readRDS("Lc_pbmc.rds")
pbmc$group_clusters <- paste0(pbmc$custom_clusters, "_", pbmc$seurat_clusters) 
colours <- c(
"#FF9A76",
"#8963BA",
"#34688D",
"#6CB9A6",
"#C1719C",
"#64966E",
"#A53D29"
)
color_plot <- colours
names(color_plot)<-unique(Idents(pbmc))
do_DimPlot(pbmc,pt.size=0.1,reduction ="umap",colors.use = color_plot) 

cell_counts <- table(Idents(pbmc))

cell_prop <- as.data.frame(cell_counts)
colnames(cell_prop) <- c("Cluster", "Cell_Number")
cell_prop$Proportion <- cell_prop$Cell_Number / sum(cell_prop$Cell_Number)

# results
cell_prop
#Cluster Cell_Number  Proportion
#1  Thrombocyte       16108 0.551002258
#2     Monocyte        3800 0.129985633
#3        VLRB+        2933 0.100328385
#4     VLRA+/C+        5356 0.183211329
#5      Basophil          122 0.004173223
#6  Granulocyte         492 0.016829719
#7 Erythorocyte         423 0.014469453