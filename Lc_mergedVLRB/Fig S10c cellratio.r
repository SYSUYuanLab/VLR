library(Seurat)
library(dplyr)
library(ggplot2)

VLRB <- readRDS("Lc_merged_VLRB.rds")
#
target_organs <- c("Ki", "PBMC")  
sub_obj <- subset(VLRB, subset = group %in% target_organs)

# 
df_prop <- sub_obj@meta.data %>%
    group_by(group, seurat_clusters) %>%
    summarise(CellCount = n(), .groups = "drop") %>%
    group_by(group) %>%
    mutate(Proportion = CellCount / sum(CellCount)) %>%
    ungroup()

# 
df_prop$seurat_clusters <- as.character(df_prop$seurat_clusters)  
df_prop$seurat_clusters <- factor(df_prop$seurat_clusters, levels = as.character(0:11))

#
organ1 <- target_organs[1]
organ2 <- target_organs[2]

df_prop <- df_prop %>%
    mutate(Proportion_plot = ifelse(group == organ1, -Proportion, Proportion))

#
ggplot(df_prop, aes(x = Proportion_plot, y = seurat_clusters, fill = group)) +
    geom_bar(stat = "identity", width = 0.7, color = "white") +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
    scale_x_continuous(labels = abs) +
    scale_fill_manual(values = c("#F08961", "#8A9CC4")) +
    theme_minimal(base_size = 14) +
    theme(
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top"
    ) +
    labs(fill = "group",
         title = paste0("Cell type composition in ", organ1, " vs ", organ2),
         x = "Proportion of cells") +
    coord_cartesian(expand = TRUE) 