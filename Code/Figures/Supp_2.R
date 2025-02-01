# load libraries
library(Seurat)
library(BuenColors)

DefaultAssay(df_integrated) <- "ADT"
Idents(df_integrated) <- "annotation"

order_cells <- rev(c( "B-Cell","CD14-Mono", "CD16-Mono", "CD4", "CD8","NK", "MAIT", "Dendritic", "dnT", "gdT", "Treg"))
levels(df_integrated) <- order_cells
df_integrated@active.ident <- factor(x = df_integrated@active.ident, levels = order_cells)

citeseq_markers <- rownames(df_integrated[["ADT"]]) ### names of the cite_seq-markers

DotPlot(df_integrated,assay = "ADT",citeseq_markers, scale=T)+RotatedAxis()+ scale_color_gradientn(colors = jdb_palette("solar_extra")) + 
  RotatedAxis()+ coord_flip() +   xlab('Markers') +  ylab('Cell Type')+ 
  coord_flip()+theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7)) +
 scale_size(breaks = c(0,25,50, 75, 100), limits=c(0,100), name = "Percentage") + labs(title="Cite-Seq Markers") +
  theme(legend.text = element_text(size = 9),  
        legend.title=element_text(size = 11),
        axis.text.x = element_text(size = 10), # Adjust x-axis tick text size
        axis.text.y = element_text(size = 10))
