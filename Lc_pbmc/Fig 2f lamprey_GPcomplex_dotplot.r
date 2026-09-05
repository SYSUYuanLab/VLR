library(dplyr)
library(ggplot2)  

## Load Seurat object
pbmc <- readRDS("lamprey_scRNA.rds")

gp <- c(
    "Leu-Gene.8815","Leu-Gene.42299","Leu-Gene.39856","Leu-Gene.25542","VLRB+-Gene.240408","Leu-Gene.50904","VLRB-3UTR+CT","Leu-Gene.112","VLRB+-Gene.225998","VLRB+-Gene.86278" 
) #Itgb3 Esm1 Gp1ba Gp1bb-3 Gp1bb-2 Gp9-22 Vlrb Csf1r Gp1bb-1 Gp9-1

colorsForDataType <- c("white","grey", "#4F2A7E")

gene <- intersect(gp, rownames(pbmc))
p <- DotPlot(pbmc, features = rev(gene))
data <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]

ggplot(data, aes(x =id , y = features.plot)) +
    geom_point(
        aes(fill = avg.exp.scaled, size = pct.exp),
        color ='black',
        shape =21,
        stroke =0.01
    ) +
    xlab("") + ylab("") +
    scale_fill_gradientn(
        colors = c("white","grey", "#4F2A7E")
    ) +
    scale_size(
        range = c(0,12),
        limits = c(0,100),
        breaks = c(0,20,40,60,80,100)
    ) +
    theme(
        text = element_text(size =20),
        panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),  
          panel.background = element_blank(),  
        axis.line = element_line(colour ="black"),
        axis.text.x=element_text(color="black",size=16,angle=9,vjust=1,hjust=1),
          axis.text.y=element_text(color="black",size=16),
        legend.position ="right",
        legend.title = element_text(size =10)
    ) +
    guides(
        size = guide_legend(
            title.position ="top",
            title.hjust =0.5,
            ncol =1,
            byrow =TRUE,
            override.aes = list(stroke =0.4)
        ),
        fill = guide_colourbar(
            title.position ="top",
            title.hjust =0.5
        )
    ) + RotatedAxis()+coord_flip()