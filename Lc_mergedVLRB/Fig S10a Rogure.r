library(ROGUE) #2.1.7
library(tidyverse)
library(Seurat)
library(ggplot2)
library(monocle)
library(tidyverse)
library(cowplot)
library(harmony)
library(ggpubr)
library(patchwork)
library(ggsci)
library(reshape2)
library(plyr)
VLRB
expr <- GetAssayData(VLRB, assay ='SCT',layer ='counts')
expr=as.matrix(expr)
meta <- VLRB@meta.data
expr <- matr.filter(expr, min.cells = 10, min.genes = 10)

ent.res <- SE_fun(expr)
head(ent.res)
SEplot(ent.res)

rogue.value <- CalculateRogue(ent.res, platform ="UMI")
rogue.value

rogue.res <- rogue(expr, labels = meta$SCT_snn_res.0.57, samples = meta$sample, platform ="UMI", span = 0.6)
rogue.res

rogue.boxplot(rogue.res)

my_colors <-   c("#346284","#5F86A0","#8AA7BA","#B5C8D1","#6BA89C","#9FC4B9",
                    "#C6D8C0","#E8CFAF","#E1A48C","#D38382","#C3A0B9","#A5A6C5")
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

rogue_mat <- rogue.res

rogue_df <- as.data.frame(rogue_mat) %>%
    tibble::rownames_to_column("sample") %>%
    pivot_longer(
        cols = -sample,
        names_to = "cluster",
        values_to = "score"
    ) %>%
    filter(!is.na(score)) %>%
    mutate(
        cluster = factor(as.numeric(cluster), levels = sort(unique(as.numeric(cluster)))),
        score = as.numeric(score)
    )
cluster_colors <- my_colors
ggplot(rogue_df, aes(x = cluster, y = score, color = cluster)) +
    geom_jitter(
        width = 0.15,
        size = 3,
        alpha = 0.9
    ) +
    stat_summary(
        fun = median,
        geom = "crossbar",
        width = 0.5,
        linewidth = 0.5,
        color = "black"
    ) +
    scale_color_manual(values = cluster_colors) +
    theme_minimal(base_size = 14) +
    labs(
        x = "VLRB subpopulation",
        y = "ROGUE Score"
    ) +
    theme(
        axis.text.x = element_text(
            hjust = 1,
            face = "bold",
            size = 12
        ),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(
            size = 14,
            face = "bold"
        ),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(
            color = "grey80",
            linewidth = 0.5
        ),
        axis.line = element_line(
            linewidth = 0.8,
            color = "black"
        )
    ) +
    ggtitle("ROGUE scores across VLRB subpopulations")