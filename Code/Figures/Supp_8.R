library(ggplot2)
library(Seurat)

Celltype_markers <- list() # Create a list to store markers for each Cell Type
for (i in 1:length(order_cells)) {
  markers_for_contrast <- list()   # Create a list to store markers for each contrast
  for (j in 1:length(patients)) {
    print(paste("Processing:", order_cells[i], "for", patients[j]))
    cell_type <- subset(x = data, subset = annotation == order_cells[i])
    Markers <- FindMarkers(cell_type, ident.1 = patients[j], ident.2 = "Healthy", min.pct = 0.5)
    Markers$Gene <- rownames(Markers)
    Markers$diffexpressed <- "Not sig."
    Markers$diffexpressed[(Markers$p_val_adj )< 0.05 ] <- "P-value"
    Markers$diffexpressed[abs(Markers$avg_log2FC )> 0.58 & Markers$p_val_adj>0.05] <- "Log(base 2) FC"
    Markers$diffexpressed[(Markers$avg_log2FC )> 0.58 & Markers$p_val_adj<0.05] <- "P-value & Log2 FC"
    Markers$diffexpressed[Markers$avg_log2FC < -0.58 & Markers$p_val_adj<0.05] <- "P-value & Log2 FC"
    Markers$plog10 <- -log10(Markers$p_val_adj)
    Markers$plog10[Markers$plog10=="Inf"] <- round(max(Markers$plog10[Markers$plog10 != Inf]))
    markers_for_contrast[[j]] <- Markers  # Store markers for this patient
  }
  Celltype_markers[[i]] <- markers_for_contrast  # Store markers for this cell type
  
}
names(Celltype_markers) <- order_cells

Volcano_plot <- function(df, Celltype) {
  ggplot(df, aes(x=avg_log2FC, y=-log10(p_val_adj), colour=diffexpressed)) +
    geom_point(size=0.5, alpha=0.5) + 
    geom_vline(xintercept = c(-0.58, 0.58), linetype="dotted", color = "black", size=0.5) +
    geom_hline(yintercept = c(log10(0.05)*(-1)), linetype="dotted", color = "black", size=0.5) +
    theme_minimal() +
    scale_color_manual(values=c(`Not sig.`="darkgrey",`P-value`="blue",`Log(base 2) FC`="darkgreen",`P-value & Log2 FC` ="red")) +
    xlim(c(-3, 3)) + 
    ylim (c(0,310))+
    ggtitle(Celltype)  +
    
    labs(x = "", y = "", color = "") +
    theme(
      legend.text =   element_text(size=10),
      axis.text.x = element_text(hjust = 1, size=11),
      axis.text.y = element_text(size=11),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10)) +
    annotate("text", x = -2, y = 220, label = paste("UP", dim(df[df$avg_log2FC > 0.58 & df$p_val_adj < 0.05, ])[1], sep = ":"), color = "black", size = 3) +
    annotate("text", x = -2, y = 190, label = paste("DOWN", dim(df[df$avg_log2FC < -0.58 & df$p_val_adj < 0.05, ])[1], sep = ":"), color = "black", size = 3)
}

#### preACLF (3) vs HC

a <- Volcano_plot(Celltype_markers$`B-Cell`[[3]], "B cells")
b <- Volcano_plot(Celltype_markers$`CD14-Mono`[[3]], "CD14 Monocytes")
c <- Volcano_plot(Celltype_markers$`CD16-Mono`[[3]], "CD16 Monocytes")
d <- Volcano_plot(Celltype_markers$CD4[[3]], "CD4+ T cells")
e <- Volcano_plot(Celltype_markers$CD8[[3]], "CD8+ T cells")
f <- Volcano_plot(Celltype_markers$NK[[3]], "NK cells")
g <- Volcano_plot(Celltype_markers$MAIT[[3]], "MAIT cells")
h <- Volcano_plot(Celltype_markers$Treg[[3]], "Tregs")

volcanos_ACLF <- ggarrange(a,b,c,d,e,f,g,h,common.legend = TRUE,  ncol = 2, nrow=4)
annotate_figure(volcanos_ACLF,top = text_grob("preACLF vs HC", face = "bold", size = 14),  left = textGrob("-Log10 Adjusted P-value", rot = 90, vjust = 1, gp = gpar(cex = 1)),
                bottom = textGrob("Fold Change", gp = gpar(cex = 1)))

