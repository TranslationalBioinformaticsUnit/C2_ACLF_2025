# load libraries
library(BuenColors)
library(Seurat)

# A) UMAP with CD14 zoom
data@meta.data$all <- as.character(data$groups)
data@meta.data$all[data@meta.data$groups=="disease"] <- "cells"
data@meta.data$all[data@meta.data$groups=="healthy"] <- "cells"

Idents(data) <- "annotation"

CD14 <- subset(x = data, idents = c("CD14-Mono"))
Idents(CD14) <- "annot2"

### UMAP with CD14 in red and rest in grey
data@meta.data$grey <- as.character(data$annotation)
data@meta.data$grey[data@meta.data$annotation=="CD14-Mono"] <- "CD14 Monocytes"
data@meta.data$grey[!data@meta.data$annotation=="CD14-Mono"] <- " "

Idents(data) <- "grey"
Umap_2colors <- DimPlot(data, raster = F, cols=c("snow3", "ivory4"), label = T) + ggtitle(NULL) + NoAxes() +
  theme( legend.position = "none") 

### Coloured CD14 umap with subclusters
Umap_CD14 <- DimPlot(CD14, group.by = "annotation", raster=F, cols = c("#d62728", "lightpink", "deeppink4")) + ggtitle(NULL) + NoAxes() +theme( legend.position = "none") 

#### Barplot for CD14s including subtypes summarise in one bar
proportions_CD14 <- data.frame(find_proportions_df(
  CD14,
  x = "all",
  fill = "annotation"
))

proportions_CD14$annotation <- as.factor(proportions_CD14$annotation)
levels(proportions_CD14$annotation) <- (c("C0" ,"C1", "C2"))

proportions_CD14 <- proportions_CD14[order(proportions_CD14$annotation),]

color_palette3 <-  (c( "#d62728", "lightpink", "deeppink4"))

stacked_barplot_CD14 <- proportions_CD14 %>%
  ggplot(aes(x = all, y = percentage_cells, group = annotation, fill = annotation)) +
  geom_col(position = position_stack(reverse = TRUE))+
  ggtitle("") + 
  coord_flip() +
  labs(x = " ", y = "Percentage of Cells (%)", fill = "") + 
  scale_fill_manual(values = color_palette3) +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    axis.text.y = element_blank(),  # Remove y-axis tick text
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )


# B) Barplot with relative proportions of CD14 by patient 

proportions_df_pat <- find_proportions_df(
  data,
  x = "patient",
  fill = "annotation"
)

CD14_prop_pat <- proportions_df_pat[proportions_df_pat$annotation=="C0" | proportions_df_pat$annotation=="C1" | proportions_df_pat$annotation=="C2", ]

barplot_CD14 <- CD14_prop_pat %>% 
  ggplot(aes(x = patient, y = percentage_cells, group = annotation, fill = annotation)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values =  c("#d62728", "lightpink", "deeppink4"))+
  labs(x = "Patient", 
       y = "Percent (%)",
       fill = "CD14+ Monocytes") +
  theme(
    # Hide panel borders and remove grid lines
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    # Change axis 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size=10),
    axis.text.y = element_text(size=10),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    legend.position = "none"
  )


# 3 C) Cleveland for CD14 subclusters   
data$orig.ident <- data@meta.data$groups

prop_test <- sc_utils(data)

proptest_res <- list()
comp <- c("SDC", "UDC","ACLF_early")
celltype <- names(table(data$annot2))

for (i in 1:length(celltype)){
  res1 <- permutation_test(
    prop_test, cluster_identity = "annot2", #predicted.id
    sample_1 = "Healthy", sample_2 = comp[i],
    sample_identity = "orig.ident")
  
  res <- c(res1@results$permutation)  # pvalue cd14
  proptest_res[[i]] <- res
}
names(proptest_res) <- comp

proptest_res_df <- data.frame(celltypes=rep(proptest_res[[1]]$clusters,3),
                              Pval= c(proptest_res[[1]]$FDR,proptest_res[[2]]$FDR,proptest_res[[3]]$FDR),
                              FC=c(proptest_res[[1]]$obs_log2FD,proptest_res[[2]]$obs_log2FD,proptest_res[[3]]$obs_log2FD),
                              Contrast= c(rep("SDC vs HC",13), rep("UDC vs HC",13), rep("preACLF vs HC",13)),
                              Color= c(rep("#FFCC00",13), rep("#FF6600",13), rep("#990000",13) ))


proptest_res_df_sub$Contrast <- factor(proptest_res_df_sub$Contrast ,levels=(c("SDC vs HC", "UDC vs HC", "preACLF vs HC")))
celltype2 <- c("C0"  ,   "C1"  ,   "C2")
proptest_res_df_sub <- proptest_res_df[proptest_res_df$celltypes %in% celltype2,]
proptest_res_df_sub$celltypes <- factor(proptest_res_df_sub$celltypes, levels=rev(c("C0"  ,   "C1"  ,   "C2")))

cleveland_Monocytes <- ggplot(proptest_res_df_sub, aes(FC, celltypes, color=Contrast)) +
  geom_point(aes(size=Pval),alpha =0.8) +
  guides(size = guide_legend(reverse=F))+
  scale_size("P.value (FDR)",range=c(4,1), breaks=c(0.0015, 0.01, 0.05, 0.1, 0.5)) +
  scale_colour_manual(values=setNames(proptest_res_df$Color,proptest_res_df$Contrast))+
  xlim(-1, 4)+ geom_vline(xintercept = c(-0.58,0.58), linetype="dotted", color = "black", size=0.5)+
  labs(x = "Log2FC", y = "CD14+ Monocytes", fill = "") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=10),
    axis.text.y = element_text(size=10),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    
    legend.title = element_text( size=12), legend.text = element_text(size = 10)
    
  )


##### D) Trajectory Analysis 

cds_short_names_df <- data.frame(gene_short_name=rownames(CD14))
rownames(cds_short_names_df) <- rownames(CD14)
cds <- new_cell_data_set(CD14@assays$RNA@data , cell_metadata =CD14@meta.data, gene_metadata=cds_short_names_df)


cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- CD14@reductions[["umap"]]@cell.embeddings
cds<- cluster_cells(cds)

cds@clusters@listData[["UMAP"]][["clusters"]] <- CD14@active.ident

cds <- learn_graph(cds)
cds <- order_cells(cds, reduction_method = "UMAP")
plot_cells(cds, color_cells_by="pseudotime", label_groups_by_cluster=F, label_branch_points=F, label_roots=F, label_leaves=F)


###### E) Dotplot 
markers_list <- c("CD14", "HLA-DRA","HLA-DRB1","LYZ","S100A8", "S100A9", "CD68","FUT4","HIF1A", "IRAK3","MERTK","TLR2", "VCAN")

CD14_Dotplot <- DotPlot(CD14, markers_list, scale=T) +
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +  RotatedAxis()+ xlab('') +  ylab('') + 
  scale_size(breaks = c(0,25,50, 75, 100), limits=c(0,100), name = "Percentage") + labs(title=" ") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 10), 
          

