library(reticulate)
use_python("~/miniconda3/envs/magic_env/bin/python3", required = TRUE) 
py_config()
library(Rmagic)
library(ggplot2)
library(readr)
library(viridis)
library(phateR)
library(corrplot)
library(dplyr)
library(purrr)
pbmc
counts_mat <-t(as.matrix(Seurat::GetAssayData(pbmc, layer = 'data')))
keep_cols <- colSums(counts_mat > 0) > 10
counts_mat <- counts_mat[,keep_cols]
ggplot() +
geom_histogram(aes(x=rowSums(counts_mat)), bins=50) +
geom_vline(xintercept = 1000, color='red') 

complex_MAGIC <- magic(counts_mat,genes = c("VLRB-3UTR+CT","VLRB+-Gene.225998","VLRB+-Gene.86278","Leu-Gene.50904","VLRB+-Gene.240408","Leu-Gene.25542","Leu-Gene.39856")) 

#before magic
ggplot(counts_mat) +
geom_point(aes(`VLRB+-Gene.225998`, `VLRB+-Gene.86278`, color=`VLRB-3UTR+CT`)) +
scale_color_viridis(option="B")

#after magic
ggplot(complex_MAGIC) +
geom_point(aes(`VLRB+-Gene.225998`, `VLRB+-Gene.86278`, color=`VLRB-3UTR+CT`)) +
scale_color_viridis(option="B")

complex_MAGIC <- magic(counts_mat,genes = c("VLRB-3UTR+CT","VLRB+-Gene.225998","VLRB+-Gene.86278","Leu-Gene.50904","VLRB+-Gene.240408","Leu-Gene.25542","Leu-Gene.39856"),t=5,init = complex_MAGIC) 

#all genes
complex_MAGIC <- magic(counts_mat, genes="all_genes",
                    knn=15,init=complex_MAGIC)

as.data.frame(complex_MAGIC)[1:5, 1:10]
complex <- complex_MAGIC$result[,c("VLRB-3UTR+CT","VLRB+-Gene.225998","VLRB+-Gene.86278","Leu-Gene.50904","VLRB+-Gene.240408","Leu-Gene.25542","Leu-Gene.39856")]

library(Hmisc)
complex_cor <- rcorr(as.matrix(complex),type = "spearman")
complex_cor

#subset
cell_group <-pbmc$custom_clusters
head(cell_group)
complex_MAGIC_gene <- complex_MAGIC$result[,c("VLRB-3UTR+CT","VLRB+-Gene.225998","VLRB+-Gene.86278","Leu-Gene.39856","Leu-Gene.50904","VLRB+-Gene.240408","Leu-Gene.25542")]
# merge
complexMAGIC_with_groups <- cbind(complex_MAGIC_gene, cell_group)

#1
complexMAGIC_with_groups$cell_group3 <- dplyr::case_when(
  complexMAGIC_with_groups$cell_group == "VLRB" ~ "VLRB",
  complexMAGIC_with_groups$cell_group == "Thrombocyte" ~ "Thrombocyte",
  TRUE ~ "Other"
)
table(complexMAGIC_with_groups$cell_group3)

#2
gene_pairs <- list(
  c("VLRB-3UTR+CT", "VLRB+-Gene.225998"), #Vlrb Gp1bb-1
  c("VLRB-3UTR+CT", "VLRB+-Gene.86278"), #Vlrb Gp9-1
  c("VLRB+-Gene.225998", "VLRB+-Gene.86278"), # Gp1bb-1 Gp9-1
  c("VLRB-3UTR+CT", "VLRB+-Gene.240408"),  #Vlrb GP1bb-2
    c("VLRB-3UTR+CT", "Leu-Gene.25542"), #Vlrb Gp1bb-3
 c("Leu-Gene.50904", "Leu-Gene.39856"), #GP9-2 Gp1ba
  c("VLRB+-Gene.240408", "Leu-Gene.39856"),  #Gp1bb-2 Gp1ba
     c("Leu-Gene.25542", "Leu-Gene.39856"),#GP1bb-3 Gp1ba
    c("Leu-Gene.50904", "VLRB+-Gene.240408"),#GP9-2 Gp1bb-2
    c("Leu-Gene.50904", "Leu-Gene.25542"),#GP9-2 Gp1bb-3
    c("VLRB+-Gene.240408", "Leu-Gene.25542"), #Gp1bb-2  Gp1bb-3
    c("VLRB-3UTR+CT", "Leu-Gene.50904") #Vlrb Gp9-2
)


plot_corr_pair <- function(pair, data){
  gene1 <- pair[1]
  gene2 <- pair[2]
  cor_results <- data %>%
    group_by(cell_group3) %>%
    summarise(
      rho = cor(.data[[gene1]],
                .data[[gene2]],
                method="spearman",
                use="complete.obs"),
      p = cor.test(.data[[gene1]],
                   .data[[gene2]],
                   method="spearman")$p.value,
      .groups="drop"
    )
  cor_results$label <- paste0(
  cor_results$cell_group3,
  ": ρ=",
  round(cor_results$rho,2),
  "\nP=",
  signif(cor_results$p,3)
    )
  cor_results <- cor_results %>%
    mutate(
      y_pos = max(data[[gene2]], na.rm=TRUE) *
        (1 - 0.1*(row_number()-1))
    )
  p <- ggplot(data,
              aes(x=.data[[gene1]],
                  y=.data[[gene2]],
                  color=cell_group3))+
    geom_point(size=0.5, alpha=0.5)+
   scale_color_manual(
  values = c(
    "VLRB" = "#E69F00",
    "Thrombocyte" = "#A6CEE3",
    "Other" = "#BDBDBD"
  )
)+
    theme_classic()+
    geom_text(
      data=cor_results,
      aes(
        x=Inf,
        y=y_pos,
        label=label,
        color=cell_group3
      ),
      hjust=1.1,
      size=3
    )+
    labs(
      title=paste(gene1,"vs",gene2),
      x=gene1,
      y=gene2,
      color="Cell group"
    )

  return(p)
}


out_dir <- "~/results"


walk(gene_pairs, function(pair){

  p <- plot_corr_pair(
    pair,
    complexMAGIC_with_groups
  )

  file_name <- paste0(
    "scatter_",
    gsub("[^A-Za-z0-9]", "_", pair[1]),
    "_vs_",
    gsub("[^A-Za-z0-9]", "_", pair[2]),
    ".pdf"
  )

  ggsave(
    filename=file.path(out_dir,file_name),
    plot=p,
    width=4,
    height=3
  )

})
write.csv(
  cor_results,
  file.path(out_dir,
            paste0(pair[1],"_vs_",pair[2],"_cor.csv")),
  row.names=FALSE
)