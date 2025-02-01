
# Load packages
library(Seurat)
library(ggVennDiagram)
library(ggvenn)
library(ComplexHeatmap)



# A) Venn Diagram 
x <- list( `CD14 vs all cell types`=comparative1$genes,
           `C2 vs C0+C1` = comparative2$genes,
           `C2 pre-ACLF vs C2 all groups` = comparative3$genes)

ggvenn(x, fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4)



# B) Heatmap

mat <- DEA_results
names(mat) <- c("pre-ACLF vs HC","pre-ACLF vs SDC", "pre-ACLF vs UDC" )

legend <- function() {
  grid.text("*: p.adj<0.05 and FC=>1.5", x = 1, y = 0.5, just = "left", gp = gpar(fontsize = 8))
  grid.text("X: FC=>1.25", x = 1, y = 0.4, just = "left", gp = gpar(fontsize = 8))
}

Heatmap(mat, name = "log2 FC",cluster_columns = F,
        row_names_gp = gpar(fontsize =6), 
        column_names_gp =  gpar(fontsize = 8),
        column_names_rot =45 ,
        row_dend_width = unit(0.6, "cm"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (sig_mat[i, j] == 1) {
            grid.text("*", x, y, gp = gpar(fontsize = 5))
          } else if (sig_mat[i, j] ==0.25 |sig_mat[i, j]==0.5 ) {
            grid.text("X", x, y, gp = gpar(fontsize = 5))
          }
        })



# C) UMAP signature with candidate genes

get_umap_signature <- function(data, color_by_pos, title){
  data <-  data[order(data[[color_by_pos]]),]
  p <- ggplot(data, aes(x=UMAP_1, y= UMAP_2, color= .data[[color_by_pos]])) + 
    geom_point(size=0.5, alpha=0.6) +
    scale_color_gradientn(colors = jdb_palette("solar_extra")) + ggtitle(title) + 
    theme_classic() + theme(legend.position = "right",panel.background = element_blank(),
                            axis.ticks = element_blank(), axis.text = element_blank())
  return(p)
}  

intersection_list <- list(candidate_genes)

tmp <- FetchData(
  AddModuleScore(data, features=intersection_list, name = 'Signature', assay='RNA'),
  vars=c('UMAP_1', 'UMAP_2', 'annotation' , 'Signature1'))

signature_all <- get_umap_signature(tmp, "Signature1", "C2 signature")

data_HC <- subset(x = CD14, subset=groups_new == c("Healthy"))
tmp_HC <- FetchData(
  AddModuleScore(data_HC, features=intersection_list, name = 'Signature', assay='RNA'),
  vars=c('UMAP_1', 'UMAP_2', 'annot2' , 'Signature1'))
signature_HC <- get_umap_signature2(tmp_HC, "Signature1", "HC")

data_SDC <- subset(x = CD14, subset=groups_new == c("SDC"))
tmp_SDC <- FetchData(
  AddModuleScore(data_SDC, features=intersection_list, name = 'Signature', assay='RNA'),
  vars=c('UMAP_1', 'UMAP_2', 'annot2' , 'Signature1'))
signature_SDC <- get_umap_signature2(tmp_SDC, "Signature1", "SDC")

data_UDC <-  subset(x = CD14, subset=groups_new == c("UDC"))
tmp_UDC <- FetchData(
  AddModuleScore(data_UDC, features=intersection_list, name = 'Signature', assay='RNA'),
  vars=c('UMAP_1', 'UMAP_2', 'annot2' , 'Signature1'))
signature_UDC <- get_umap_signature2(tmp_UDC, "Signature1", "UDC")

data_preACLF <-  subset(x = CD14, subset=groups_new == c("ACLF"))
tmp_preACLF <- FetchData(
  AddModuleScore(data_preACLF, features=intersection_list, name = 'Signature', assay='RNA'),
  vars=c('UMAP_1', 'UMAP_2', 'annot2' , 'Signature1'))
signature_preACLF <- get_umap_signature2(tmp_preACLF, "Signature1", "pre-ACLF")

signature_groups <- ggarrange(signature_HC,signature_SDC, signature_UDC, signature_preACLF)



# D) Boxplot

ggplot(data = geneSet_df, aes(x=Group, y=GeneSet)) + geom_boxplot(aes(fill=Group))+
  scale_fill_manual(values=c("darkgreen", "gold1", "darkorange1", "firebrick4")) + 
  labs(x = "", y = "", title = "Bulk RNA-seq data (N=19) ", fill="") +
       theme(legend.position="none",
             legend.margin = margin(t = -15, r = 0, b = 0, l = 0),
             axis.text=element_text(size=10),
             axis.text.x=element_text(size=10, angle=45,hjust = 1),
             axis.text.y=element_text(size=10),plot.title = element_text(face = "bold", hjust = 0.5, size = 11))


