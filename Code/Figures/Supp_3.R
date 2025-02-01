# load libraries
library(Seurat)
library(BuenColors)
library(scProportionTest)

# prepare object
DefaultAssay(df_integrated) <- "ADT"
Idents(df_integrated) <- "annotation"
order_cells <- rev(c( "B-Cell","CD14-Mono", "CD16-Mono", "CD4", "CD8","NK", "MAIT", "Dendritic", "dnT", "gdT", "Treg"))
levels(df_integrated) <- order_cells
df_integrated@active.ident <- factor(x = df_integrated@active.ident, levels = order_cells)

# A) Annotated umap
Umap_data <- DimPlot(df_integrated, group.by = "annotation", raster=F, label = T) + ggtitle("CITE-seq") + NoAxes() + 
  theme(legend.position = "bottom", legend.text =   element_text(size=8), legend.key.size = unit(0.5, "cm"))

# B) Dotplot
citeseq_markers <- rownames(df_integrated[["ADT"]]) ### names of the cite_seq-markers

DotPlot(df_integrated,assay = "ADT",citeseq_markers, scale=T)+RotatedAxis()+ scale_color_gradientn(colors = jdb_palette("solar_extra")) + 
  RotatedAxis()+ coord_flip() +   xlab('Markers') +  ylab('Cell Type')+ 
  coord_flip()+theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7)) +
 scale_size(breaks = c(0,25,50, 75, 100), limits=c(0,100), name = "Percentage") + labs(title="Cite-Seq Markers") +
  theme(legend.text = element_text(size = 9),  
        legend.title=element_text(size = 11),
        axis.text.x = element_text(size = 10), # Adjust x-axis tick text size
        axis.text.y = element_text(size = 10))

# C) Barplot with proportions
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
  df_integrated,
  x = "groups",
  fill = "annotation"
)
stacked_barplot_cells <- plot_stacked_barplot(
  df_integrated,
  x = "groups",
  fill = "annotation",
  colors = color_palette
)


# D) Cleveland plot 
prop_test <- sc_utils(df_integrated)


proptest_res <- list()
comp <- c("SDC", "UDC","pre_aclf")
celltype <- names(table(df_integrated$annotation))

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

