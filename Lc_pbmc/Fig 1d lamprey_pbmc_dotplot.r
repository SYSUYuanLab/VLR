library(Seurat)
library(ggplot2) 

## 1. Load Seurat object
pbmc <- readRDS("lamprey_scRNA.rds")

# 2. Marker genes for DotPlot
marker_genes <- list(
  Thrombocyte = c("Leu-Gene.8815","Leu-Gene.42299","Leu-Gene.65742"),  #Itgb3 Esm1 Pxn
  Monocyte = c("Leu-Gene.17619","Leu-Gene.59328","Leu-Gene.16442","Leu-Gene.15461"),   #Mmp20/13  Mmp2 Mmp19 Ltc4s
  VLRB = c("VLRB-3UTR+CT","Leu-Gene.112","Leu-Gene.45313"),   #Vlrb Csf1r Pnn
  `T-like cell` = c("Leu-Gene.19668","Leu-Gene.21326","Leu-Gene.80706"), #Vlrc Vlra Apex1
  Basophil = c("Leu-Gene.60196","Leu-Gene.34442","Leu-Gene.70721"),  # CD22 SULT2B1 CD63
  Granulocyte = c("Leu-Gene.51667","VLRB+-Gene.254297","Leu-Gene.29730"),   # C1qtnf9 Lrp1 Aif1
  Erythrocyte = c("Leu-Gene.49859","Leu-Gene.37979","Leu-Gene.87153")   #GLOB1 Cygb Epb41
)

# 3. DotPlot
DotPlot(pbmc,features = marker_genes,assay = "SCT",scale = TRUE
) +
  scale_color_gradientn(colors = c("white", "grey", "#4F2A7E")) +
  scale_size(range = c(0,10),limits = c(0,100)) +
  theme_bw() +
  theme(panel.grid = element_blank(),panel.background = element_blank(),legend.position = "right",strip.text = element_text(size = 12,color = "black"),axis.text.x = element_text(color = "black",size = 12,angle = 45,vjust = 1,hjust = 1),
axis.text.y = element_text(color = "black",size = 12),
    legend.text = element_text(size = 10,color = "black"),
    legend.title = element_text(size = 10,color = "black" )
  )