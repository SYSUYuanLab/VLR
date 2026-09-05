VLRB <- readRDS("Lc_merged_VLRB.rds")
VLRB
library(Seurat)
library(ggplot2)
library(dplyr)
my_colors <-   c("#346284","#5F86A0","#8AA7BA","#B5C8D1","#6BA89C","#9FC4B9",
                    "#C6D8C0","#E8CFAF","#E1A48C","#D38382","#C3A0B9","#A5A6C5")

DimPlot(VLRB,reduction = "umap",pt.size = 1,label = TRUE,group.by = "group",cols = c("Set2")) #Fig.6A
DimPlot(VLRB,reduction = "umap",pt.size = 1,label = TRUE,split.by = "group",cols = my_colors) #Fig.6B_top
table(VLRB$group)
#Gi   Ki PBMC 
 #71  971 3056 
 
 #Fig.6b
library(ggplot2)
library(dplyr)

count_df <- VLRB@meta.data %>%
  filter(group %in% c("Ki", "PBMC")) %>%
  group_by(group, seurat_clusters) %>%
  summarise(count = n(), .groups = "drop")

count_df <- count_df %>%
  mutate(count_plot = ifelse(group == "PBMC", -count, count))

count_df$seurat_clusters <- factor(count_df$seurat_clusters, 
                                   levels = unique(count_df$seurat_clusters))


ggplot(count_df, aes(x = seurat_clusters, y = count_plot, fill = seurat_clusters)) +
  geom_bar(stat="identity", width=0.7) +
  scale_y_continuous(labels = abs, 
                     limits = c(-max(count_df$count)*1.1, max(count_df$count)*1.1)) +
  scale_fill_manual(values = my_colors) +
  theme_minimal(base_size = 14) +
  labs(x = "B cell subcluster", y = "Cell count", fill = "Tissue") +
  theme(
    axis.text.x = element_text(size=12, angle=45, hjust=1),
    axis.text.y = element_text(size=12),
    axis.title = element_text(size=13),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust=0.5, face="bold")
  ) +
  ggtitle("B cell subcluster distribution: Kidney vs PBMC")

  scale_fill_manual(values = c("#346284","#5F86A0","#8AA7BA","#B5C8D1","#6BA89C","#9FC4B9",
                    "#C6D8C0","#E8CFAF","#E1A48C","#D38382","#C3A0B9","#A5A6C5"))