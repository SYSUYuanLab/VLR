library(Biobase)
library(monocle) #dplyr 1.1.0
library(Seurat)
library(dplyr)
VLRB <- readRDS("Lc_merged_VLRB.rds")
VLRB <- subset(VLRB, subset = group != "Gi")

DefaultAssay(VLRB) <- "RNA"
VLRB <- JoinLayers(VLRB)
expr_matrix <- GetAssayData(VLRB, assay = "RNA", layer = "counts")
sample_sheet <- VLRB@meta.data 
gene_annotation=data.frame(gene_short_name=rownames(VLRB))
rownames(gene_annotation) <- rownames(VLRB)
pd<-new("AnnotatedDataFrame",data=sample_sheet)
fd<-new("AnnotatedDataFrame",data=gene_annotation)
cds<-newCellDataSet(expr_matrix,phenoData=pd,featureData=fd,expressionFamily=negbinomial.size())
    
cds
cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),
num_cells_expressed >= 10))
#1
DefaultAssay(VLRB) <- "SCT"
clusters <- unique(VLRB$seurat_clusters)
VLRB <- PrepSCTFindMarkers(VLRB)
top_genes_list <- list()
for (cl in clusters) {
  markers <- FindMarkers(
    object = VLRB,
    ident.1 = cl, 
    logfc.threshold = 0.25
  )
  markers <- markers %>% filter(p_val_adj <= 0.05)
  markers <- markers %>% arrange(p_val_adj)

  n_top <- min(100, nrow(markers))
  top_genes <- rownames(markers)[1:n_top]
  top_genes_list[[as.character(cl)]] <- top_genes
}

ordering_genes <- unique(unlist(top_genes_list)) 
length(ordering_genes)
cds<-setOrderingFilter(cds,ordering_genes)
plot_ordering_genes(cds)

    
#2
plot_pc_variance_explained(cds) 
cds<-reduceDimension(cds,
                     max_components=2,
                     num_dim=8,
                     reduction_method='DDRTree',
                     residualModelFormulaStr="~sample+group",
                     verbose=F)
    
#3
cds<-orderCells(cds)
plot_cell_trajectory(cds, color_by = "State")+facet_wrap(~State,nrow=1)

##4
cds<-orderCells(cds,root_state=5)

plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE)

p <-plot_cell_trajectory(cds, color_by = "seurat_clusters")+facet_wrap(~seurat_clusters,nrow=1)
p+scale_color_manual(values =   c("0"="#346284","1"="#5F86A0","2"="#8AA7BA","3"="#B5C8D1","4"="#6BA89C","5"="#9FC4B9",     "6"="#C6D8C0","7"="#E8CFAF","8"="#E1A48C","9"="#D38382","10"="#C3A0B9","11"="#A5A6C5"))

ggplot(pData(cds), aes(Pseudotime, colour = seurat_clusters,fill = seurat_clusters) )+geom_density(bw =0.5, size =1, alpha =0.5) +theme_classic()
library(ggplot2)
pseudotime_df <- data.frame(
    pseudotime = pData(cds)$Pseudotime,
    cluster = pData(cds)$seurat_clusters
)
ggplot(pseudotime_df, aes(x = cluster, y = pseudotime, fill = cluster)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.8, color = NA) +
    geom_jitter(width = 0.25, size = 0.6, alpha = 0.5) +
    labs(
        x = "Seurat Cluster",
        y = "Pseudotime",
        title = "Distribution of pseudotime across Seurat clusters"
    ) +scale_fill_manual(values = c("0"="#346284","1"="#5F86A0","2"="#8AA7BA","3"="#B5C8D1","4"="#6BA89C","5"="#9FC4B9",     "6"="#C6D8C0","7"="#E8CFAF","8"="#E1A48C","9"="#D38382","10"="#C3A0B9","11"="#A5A6C5")) +
    theme_classic(base_size = 14) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

VLRB@meta.data$Pseudotime <- cds@phenoData@data$Pseudotime
mydata<- FetchData(VLRB,vars = c("umap_1","umap_2","Pseudotime"))
p <- ggplot(mydata,aes(x = umap_1,y =umap_2,colour = Pseudotime))+
    geom_point(size = 1)
p4 <- p + theme_bw() + theme(panel.border = element_blank(), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), 
                             axis.line = element_line(colour = "black"))

#5 BEAM one method
BEAM_res <- BEAM(cds, branch_point = 1, progenitor_method = "duplicate") #or expression or cluster y
table(pData(cds[, row.names(BEAM_res)])$State)

cds_subset <- buildBranchCellDataSet(cds, 
                                     branch_point = 1, 
                                     progenitor_method = "duplicate")
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
cds_heatmap <- plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,qval < 0.0001)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T,
                            return_heatmap=T) 
     



#6 BEAM gene
genes <- row.names(subset(fData(cds),
          gene_short_name %in% c("CDA2","VLRB-3UTR+CT","VLRB+-Gene.225998","VLRB+-Gene.86278")))
plot_genes_branched_pseudotime(cds[genes,],
                       branch_point = 1,
                       color_by = "Pseudotime",
                       ncol = 1)+scale_color_viridis_c(option = "C")
                                     
#7 single gene
library(ggsci)
pData(cds)$CDA2 = log2( exprs(cds)['CDA2',]+1)
plot_cell_trajectory(cds, color_by = "CDA2",cell_size=0.5)  + scale_color_gsea() #600 600 2026.01.28_Fig5_Gene.86278BEAM

pData(cds)$VLRB = log2( exprs(cds)['VLRB-3UTR+CT',]+1)
plot_cell_trajectory(cds, color_by = "VLRB",cell_size=0.5)  + scale_color_gsea()

pData(cds)$Gene.225998 = log2( exprs(cds)['VLRB+-Gene.225998',]+1)
plot_cell_trajectory(cds, color_by = "Gene.225998",cell_size=0.5)  + scale_color_gsea()

pData(cds)$Gene.86278 = log2( exprs(cds)['VLRB+-Gene.86278',]+1)
plot_cell_trajectory(cds, color_by = "Gene.86278",cell_size=0.5)  + scale_color_gsea()