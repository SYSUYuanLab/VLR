library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidydr)
library(stringr)
library(viridis)
library(scCustomize)
gene <- c("XLOC-014698","XLOC-018679","XLOC-006008","XLOC-002892","XLOC-010681","XLOC-015838","XLOC-004768","XLOC-014940","XLOC-009447") #EbVLRA EbVLRC EbVLRB GP1BB-1 GP1BB-2 GP9-1 EbGP1BA GP1BB-3 GP9-2 
i=1
plots=list()
#featureplot
for (i in 1:length(gene)){
    plots[[i]]=FeaturePlot_scCustom(pbmc,
                                    colors_use=colorRampPalette(c("#3288BD","white","#D53E4F"))(50),
                                    features=gene[i],label=FALSE)+NoAxes()#+
    #theme(panel.border=element_rect(fill=NA,color="black",
    #size=1.5,linetype="solid"))
}

library(patchwork)
p<-wrap_plots(plots,ncol=3);
p 


#dotplot
gene <- intersect(gene, rownames(pbmc))
p <- DotPlot(pbmc, features = rev(gene))
data <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
table(data$id)

ggplot(data, aes(x =id , y = features.plot)) +
    geom_point(
        aes(fill = avg.exp.scaled, size = pct.exp),
        color ='black',
        shape =21,
        stroke =0.01
    ) +
    xlab("") + ylab("") +
    scale_fill_gradientn(
        colors=c("#5bacdb","#62b8e3","#c5e1f0","white","#f0b4c9","#ea87aa","#ea7ea3")
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