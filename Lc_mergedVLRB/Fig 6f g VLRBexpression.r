library(Seurat)
VLRB <- subset(VLRB, subset = group != "Gi")
avg_exp <- AverageExpression(
  VLRB,
  features = c("VLRB-3UTR+CT","VLRB+-Gene.225998","VLRB+-Gene.86278"), #Vlrb Gp1bb-1 Gp9-1
  assays = "SCT",
  slot = "data"
)$SCT

avg_exp_scaled <- t(scale(t(avg_exp)))

library(pheatmap)

library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(
  c(-1, 0, 1),
  c("#4575B4", "white", "#D73027")
)

Heatmap(
  mat,
  col = col_fun,
  name = "Spearman",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9)
)


#addmodulescore
VLRB_genes <- list(c("VLRB-3UTR+CT","VLRB+-Gene.225998","VLRB+-Gene.86278"))
VLRB <- AddModuleScore(VLRB, features = VLRB_genes, name = "VLRB_score")

library(Seurat)
library(dplyr)

df <- data.frame(
  BCR_score = VLRB@meta.data$VLRB_score1,
  cluster = VLRB@meta.data$seurat_clusters
)

clusters <- unique(df$cluster)

result <- data.frame(
  cluster = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (cl in clusters) {
  x <- df$BCR_score[df$cluster == cl]
  y <- df$BCR_score[df$cluster != cl]
  
  test <- wilcox.test(x, y, alternative = "greater")
  
  result <- rbind(result, data.frame(cluster = cl, p_value = test$p.value))
}

result <- result %>% mutate(significant = ifelse(p_value < 0.0001, TRUE, FALSE))

result
library(ggplot2)
library(dplyr)

df <- df %>%
  mutate(cluster = reorder(cluster, BCR_score, FUN = median)) 

my_colors <-   c("0"="#346284","1"="#5F86A0","2"="#8AA7BA","3"="#B5C8D1","4"="#6BA89C","5"="#9FC4B9",     "6"="#C6D8C0","7"="#E8CFAF","8"="#E1A48C","9"="#D38382","10"="#C3A0B9","11"="#A5A6C5")
library(ggplot2)

#boxplot
ggplot(df, aes(x = cluster, y = BCR_score, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  labs(x = "B cell subtype", y = "BCR Module Score") + scale_fill_manual(values = my_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 1))

library(dplyr)
cluster_median <- df %>%
  group_by(cluster) %>%
  summarise(median_score = median(BCR_score, na.rm = TRUE)) %>%
  ungroup()