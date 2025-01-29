library(Seurat)
library(ggplot2)
library(dplyr)
library(scProportionTest)
library(grid)

data <- readRDS('./00.Data/seurat_integrated.rds')

##### A) UMAP annotated with all cells ##### 
color_palette <-  c("#1f77b4", "#d62728","#2ca02c",  "#ff7f0e", "#9467bd", "#636363","#bcbd22","#17becf", "#e377c2" ,"#8c564b","#ad494a")

Idents(data) <- "annotation"

order_cells <- c( "B-Cell","CD14-Mono", "CD16-Mono", "CD4", "CD8","NK", "MAIT", "Dendritic", "dnT", "gdT", "Treg")
levels(data) <- order_cells
data@active.ident <- factor(x = data@active.ident, levels = order_cells)

Umap_data <- DimPlot(data, group.by = "annotation", raster=F, cols = color_palette, label = T) + ggtitle(NULL) + NoAxes() + 
  theme(legend.text =   element_text(size=8), legend.key.size = unit(0.5, "cm"))

##### B) Bar with proportions of cell types by group ##### 
find_proportions_df <- function(seurat_obj, x, fill) {
  df <- seurat_obj@meta.data %>%
    dplyr::select(x, fill) %>%
    group_by(.data[[x]], .data[[fill]]) %>%
    summarise(n_cells = n()) %>%
    ungroup() %>%
    group_by(.data[[x]]) %>%
    mutate(n_cells_total = sum(n_cells)) %>%
    ungroup() %>%
    mutate(percentage_cells = round(n_cells / n_cells_total * 100, 3))
  df
}

plot_stacked_barplot <- function(df, x, fill, colors) {
  p <- df %>%
    ggplot(aes_string(x, "percentage_cells", fill = fill)) +
    geom_bar(stat="identity")+
    ggtitle("Cell type proportion") +
    labs(x = "Group", y = "Cells (%)", fill = "") + 
    scale_fill_manual(values = colors) +
    theme(
      axis.title.y = element_text(size = 10),
      axis.text=element_text(size=10),
      axis.text.x=element_text(angle=45, hjust=1,size=8),
      axis.title.x = element_text(vjust = 2, size=10),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size=11),
    ) 
  p
}

proportions_df_cells <- find_proportions_df(
  data,
  x = "groups_new",
  fill = "annotation"
)

str(proportions_df_cells)
levels(proportions_df_cells$groups_new) <- list(HC="Healthy",SDC= "SDC", UDC="UDC",preACLF="ACLF_early")
proportions_df_cells$groups_new <- factor(proportions_df_cells$groups_new, levels= (c("HC"    ,  "SDC"  ,   "UDC"  ,   "preACLF")), ordered = T)
levels(proportions_df_cells$annotation) <- list(`B cells`="B-Cell", `Classical monocytes`="CD14-Mono", `Non-classical monocytes`="CD16-Mono" ,
                                                `CD4+ T cells`= "CD4" ,`CD8+ T cells`="CD8",`NK cells`= "NK",`MAIT cells`="MAIT" ,    
                                                `Dendritic cells`="Dendritic" , `dnT cells`= "dnT", `gdT cells`="gdT", `Regulatory T cells`="Treg"   )


stacked_barplot_cells <- plot_stacked_barplot(
  proportions_df_cells,
  x = "groups_new",
  fill = "annotation",
  colors = color_palette
)


##### C) Cleveland to compare cell proportions between AD and HC ##### 
data$orig.ident <- data@meta.data$groups_new
table(data$orig.ident)
(table(data$annotation))

prop_test <- sc_utils(data)


proptest_res <- list()
comp <- c("SDC", "UDC","ACLF_early")
celltype <- names(table(data$annotation))

for (i in 1:length(celltype)){
  res1 <- permutation_test(
    prop_test, cluster_identity = "annotation", 
    sample_1 = "Healthy", sample_2 = comp[i],
    sample_identity = "orig.ident")
  
  res <- c(res1@results$permutation)  
  proptest_res[[i]] <- res
}
names(proptest_res) <- comp

proptest_res_df <- data.frame(celltypes=rep(proptest_res[[1]]$clusters,3),
                              Pval= c(proptest_res[[1]]$FDR,proptest_res[[2]]$FDR,proptest_res[[3]]$FDR),
                              FC=c(proptest_res[[1]]$obs_log2FD,proptest_res[[2]]$obs_log2FD,proptest_res[[3]]$obs_log2FD),
                              Contrast= c(rep("SDC vs HC",11), rep("UDC vs HC",11), rep("preACLF vs HC",11)),
                              Color= c(rep("#FFCC00",11), rep("#FF6600",11), rep("#990000",11) ))

table(data$annotation)/dim(df_integrated)[2]*100 # We keep those cells that represent >1 % of the global population:
celltype2 <- c("B-Cell", "CD14-Mono", "CD16-Mono", "CD4" , "CD8", "MAIT", "Treg", "NK" )
proptest_res_df_sub <- proptest_res_df[proptest_res_df$celltypes %in% celltype2,]
proptest_res_df_sub$celltypes <- factor(proptest_res_df_sub$celltypes, levels=rev(c("B-Cell", "CD14-Mono", "CD16-Mono", "CD4" , "CD8", "NK" , "MAIT", "Treg")))

proptest_res_df_sub$Contrast <- factor(proptest_res_df_sub$Contrast ,levels=(c("SDC vs HC", "UDC vs HC", "preACLF vs HC")))
cleveland <- ggplot(proptest_res_df_sub, aes(FC, celltypes, color=Contrast)) +
  ggtitle("Differences in cell type proportion") +
  geom_point(aes(size=Pval),alpha =0.8) +
  guides(size = guide_legend(reverse=F))+
  scale_size("P.value (FDR)",range=c(4,1), breaks=c(0.0015, 0.01, 0.05, 0.1, 0.5)) +
  scale_colour_manual(values=setNames(proptest_res_df$Color,proptest_res_df$Contrast))+
  xlim(-4, 4)+ geom_vline(xintercept = c(-0.58,0.58), linetype="dotted", color = "black", size=0.5)+
  labs(x = "Log2FC", y = " ", fill = "") + 
  theme(
    axis.title =element_text(size=10),
    axis.text.x=element_text(angle=45, hjust=1,size=8),
    axis.text.y=element_text(size=8),
    plot.title = element_text(face = "bold", size=11),
    legend.text =   element_text(size=8), legend.key.size = unit(0.4, "cm"), legend.title = element_text(size=9))


##### D) Number of DEA per cell type #####
Idents(data) <- "groups_new"

order_cells <- c( "B-Cell","CD14-Mono", "CD16-Mono", "CD4", "CD8","NK", "MAIT", "Dendritic", "dnT", "gdT", "Treg")
patients <- c("Healthy","SDC","UDC")

Celltype_markers <- list() # Create a list to store markers for each Cell Type
for (i in 1:length(order_cells)) {
  markers_for_contrast <- list()   # Create a list to store markers for each contrast
  for (j in 1:length(patients)) {
    print(paste("Processing:", order_cells[i], "for", patients[j]))
    cell_type <- subset(x = data, subset = annotation == order_cells[i])
    Markers <- FindMarkers(cell_type, ident.1 = "ACLF_early", ident.2 = patients[j], min.pct = 0.5)
    Markers$Gene <- rownames(Markers)
    #Markers <- subset(Markers, !startsWith(Gene, "RP"))
    Markers$diffexpressed <- "Not sig."
    Markers$diffexpressed[(Markers$p_val_adj )< 0.05 ] <- "Adj.P-value"
    Markers$diffexpressed[abs(Markers$avg_log2FC )>= 0.58 & Markers$p_val_adj>0.05] <- "Log2FC"
    Markers$diffexpressed[(Markers$avg_log2FC )>= 0.58 & Markers$p_val_adj<0.05] <- "Adj.P-value & Log2FC +"
    Markers$diffexpressed[Markers$avg_log2FC <= -0.58 & Markers$p_val_adj<0.05] <- "Adj.P-value & Log2FC -"
    Markers$plog10 <- -log10(Markers$p_val_adj)
    Markers$plog10[Markers$plog10=="Inf"] <- round(max(Markers$plog10[Markers$plog10 != Inf]))
    markers_for_contrast[[j]] <- Markers  # Store markers for this patient
  }
  Celltype_markers[[i]] <- markers_for_contrast  # Store markers for this cell type
}
names(Celltype_markers) <- order_cells

DEA_numbers <- list()
for(i in 1:length(order_cells)){
  Up <- unique(c(Celltype_markers[[i]][[1]]$Gene[Celltype_markers[[i]][[1]]$diffexpressed=="Adj.P-value & Log2FC +"],
                 Celltype_markers[[i]][[2]]$Gene[Celltype_markers[[i]][[2]]$diffexpressed=="Adj.P-value & Log2FC +"],
                 Celltype_markers[[i]][[3]]$Gene[Celltype_markers[[i]][[3]]$diffexpressed=="Adj.P-value & Log2FC +"]))
  Down <- unique(c(Celltype_markers[[1]][[1]]$Gene[Celltype_markers[[i]][[1]]$diffexpressed=="P-value & Log2FC -"],
                   Celltype_markers[[i]][[2]]$Gene[Celltype_markers[[i]][[2]]$diffexpressed=="Adj.P-value & Log2FC -"],
                   Celltype_markers[[i]][[3]]$Gene[Celltype_markers[[i]][[3]]$diffexpressed=="Adj.P-value & Log2FC -"]))
  
  DEA <- c(length(Up), length(Down))
  DEA_numbers[[i]] <- DEA
}
names(DEA_numbers) <- order_cells
DEA_numbers <- DEA_numbers[-(8:10)]
DEA_numbers <- do.call(rbind.data.frame, DEA_numbers)
rownames(DEA_numbers) <- order_cells[-(8:10)]
names(DEA_numbers) <- c("Up", "Down")
DEA_numbers <- stack(DEA_numbers)
DEA_numbers$`Cell Type` <- rep(order_cells[-(8:10)], 2)

DEA_numbers_plot <- ggplot(DEA_numbers,aes(`Cell Type`,values, fill = ind)) +
  geom_bar(stat="identity")+
  ggtitle("Number of DEGs") +
  labs(x = "Cell type", y = "Gene count", fill = "") + 
  scale_fill_manual(values = c("Up" = "red2", "Down" = "darkblue")) +  # Specify colors for ind
  theme(legend.position="top",
        legend.justification="right",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,5,-5,-5),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.4, "cm"),  
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10,vjust=3),
        axis.text=element_text(size=10),
        axis.text.x=element_text(angle=45, hjust=1, size=8),
        axis.line = element_line(color = "black", size = 0.5),
        plot.title = element_text(face = "bold", size=11))


##### E) Volcano Plots between preACLF and HC in all cell types ##### 
Volcano_plot <- function(df, Celltype) {
  ggplot(df, aes(x=avg_log2FC, y=-log10(p_val_adj), colour=diffexpressed)) +
    geom_point(size=0.5, alpha=0.5) + 
    guides(color = guide_legend(override.aes = list(size = 2)))+
    geom_vline(xintercept = c(-0.58, 0.58), linetype="dotted", color = "black", size=0.5) +
    geom_hline(yintercept = c(log10(0.05)*(-1)), linetype="dotted", color = "black", size=0.5) +
    theme_minimal() +
    scale_color_manual(values=c(`Not sig.`="darkgrey",`Adj.P-value`="chartreuse4",`Log2FC`="darkorange",`Adj.P-value & Log2FC +` ="darkblue",`Adj.P-value & Log2FC -` ="red2")) +
    xlim(c(-3, 3)) + 
    ylim (c(0,310))+
    ggtitle(Celltype)  +
    
    labs(x = "", y = "", color = "") +
    theme(
      legend.text =   element_text(size=8),
      legend.key.size = unit(0.3, "cm"),
      axis.text.x = element_text(hjust = 2, size=8),
      axis.text.y = element_text(size=8),
      axis.title.y = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(size = 10, face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 8)) 
    #annotate("text", x = -2, y = 220, label = paste("UP", dim(df[df$avg_log2FC > 0.58 & df$p_val_adj < 0.05, ])[1], sep = ":"), color = "black", size = 2) +
    #annotate("text", x = -2, y = 190, label = paste("DOWN", dim(df[df$avg_log2FC < -0.58 & df$p_val_adj < 0.05, ])[1], sep = ":"), color = "black", size = 2)
}


#### ACLF vs UDC
a <- Volcano_plot(Celltype_markers$`B-Cell`[[3]], "B cells")
b <- Volcano_plot(Celltype_markers$`CD14-Mono`[[3]], "CD14 Monocytes")
c <- Volcano_plot(Celltype_markers$`CD16-Mono`[[3]], "CD16 Monocytes")
d <- Volcano_plot(Celltype_markers$CD4[[3]], "CD4+ T cells")
e <- Volcano_plot(Celltype_markers$CD8[[3]], "CD8+ T cells")
f <- Volcano_plot(Celltype_markers$NK[[3]], "NK cells")
g <- Volcano_plot(Celltype_markers$MAIT[[3]], "MAIT cells")
h <- Volcano_plot(Celltype_markers$Treg[[3]], "Tregs")

volcanos_ACLF_UDC <- ggarrange(a,b,c,d,e,f,g,h,common.legend = TRUE,  ncol = 4, nrow=2)
volcanos_ACLF_UDC <- annotate_figure(volcanos_ACLF_UDC, top = text_grob("pre-ACLF vs UDC", face = "bold", size = 12),  
                           left = textGrob("-Log10 Adj.P-value", rot = 90, vjust = 2),
                            bottom = textGrob("Log2FC",vjust = -1.5, hjust = 0), fig.lab.size = 8)


pdf("Volcanos_alcf_UDC.pdf", height=4)
volcanos_ACLF_UDC <- ggarrange(a,b,c,d,e,f,g,h,common.legend = TRUE,  ncol = 4, nrow=2, legend = "bottom")
volcanos_ACLF_UDC<- annotate_figure(volcanos_ACLF_UDC,top = text_grob("pre-ACLF vs UDC", face = "bold", size = 12),  left = textGrob("-Log10 Adjusted P-value", rot = 90, vjust = 1, gp = gpar(cex = 1)),
                bottom = textGrob("Fold Change", gp = gpar(cex = 1)))
dev.off()
