library(SCpubr)
library(ggplot2)

VLRAC <-  c("XLOC-014698","XLOC-018679") # VLRA VLRC
VLRB <- c("XLOC-006008","XLOC-001459") #VLRB CD9
Thrombocyte <- c("XLOC-004768","XLOC-013491","XLOC-006849") #GP1BA SDPR TIMP4
Erythrocyte <- c("alas1","XLOC-019661") #alas1 GATA2
Monocyte <- c("mmp2","XLOC-011260","XLOC-005219") #IRF8,MMP19
genes_to_check=list(
    `VLRA+/C+` =VLRAC, 
    Thrombo.=Thrombocyte,
    `VLRB+`=VLRB,
    Mono. =Monocyte,
    Erythro.=Erythrocyte
)
library(stringr)
DotPlot(pbmc,
    features=genes_to_check,
    scale=T,assay='SCT')+
    theme_bw()+
 #   scale_color_continuous(low="grey",high="darkred")+ 
 scale_color_gradientn(colors = c("white","grey", "#4F2A7E")) +
    theme(
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          legend.position="right",legend.box="vertical",
          legend.margin=margin(t=0,unit='cm'),
          strip.text = element_text(size = 12, color = "black"),
          axis.text.x=element_text(color="black",size=12,angle=45,vjust=1,hjust=1),
          axis.text.y=element_text(color="black",size=12),
          legend.text=element_text(size=10,color="black"),
          legend.title=element_text(size=10,color="black"))+scale_size(range=c(0,10))  