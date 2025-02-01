# load libraries
library(Seurat)

Markers_C2 <- FindMarkers(CD14, ident.1 = "C2",  min.pct = 0.75,logfc.threshold=0.58)

Markers_order <- Markers_C2[order(Markers_C2$avg_log2FC),]

top_markers <- head(rownames(Markers_order),50)
tail_markers <-  tail(rownames(Markers_order),50)

markers <- c(top_markers,tail_markers)
CD14 <- ScaleData(object = CD14, features = c(top_markers,tail_markers))

CD14_heatmap <- DoHeatmap(CD14, features = c(top_markers,tail_markers), label = T , group.colors= c("#d62728", "lightpink", "deeppink4")) + 
   scale_fill_gradient2( low = (c("#3361A5", "#248AF3", "#14B3FF", "#88CEEF")), mid = "#C1D5DC", high =(c("#EAD397", "#FDB31A" ,"#E42A2A" ,"#A31D1D")), midpoint = 0, guide = "colourbar", aesthetics = "fill")



