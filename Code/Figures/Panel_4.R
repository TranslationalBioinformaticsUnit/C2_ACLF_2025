# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(enrichplot)
library(clusterProfiler)

Idents(CD14) <- "annotation"

Markers <- FindAllMarkers(CD14, ident.1 = "C2", min.pct = 0.5)
Markers$genes <- rownames(Markers)

genes_all <- Markers$avg_log2FC
names(genes_all) <- Markers$genes
genes_all = sort(genes_all, decreasing = TRUE)
head(genes_all)

go_GSEA <- gseGO(genes_all,
                 OrgDb = "org.Hs.eg.db", 
                 keyType = 'SYMBOL',
                 ont = "BP",
                 minGSSize = 5, maxGSSize = 250, pvalueCutoff = 1)

# A) Running score and preranked list of GSEA result
go_GSEA@result$Description <- paste(go_GSEA@result$ID, go_GSEA2@result$Description, sep = " | ")

a <- print(gseaplot2(go_GSEA, geneSetID = 1,  base_size = 8, title = go_GSEA@result$Description[1], rel_heights = c(0.25, 0.1, 0.1)))
b <- print(gseaplot2(go_GSEA, geneSetID = 2,  base_size = 8, title = go_GSEA@result$Description[2], rel_heights = c(0.25, 0.1, 0.1)))
c <- print(gseaplot2(go_GSEA, geneSetID = 16,  base_size = 8, title = go_GSEA@result$Description[3], rel_heights = c(0.1, 0.1, 0.1)))


gseas <- ggarrange(a,b,c, common.legend = TRUE, ncol = 1)


# B) Boxplots Cell exhaustion per cluster
gsvaplot_cluster <- function(Celltype, features_list, title){
  tmp <- FetchData(
    AddModuleScore(Celltype, features=features_list, name = 'Signature', assay='RNA'),
    vars=c('UMAP_1', 'UMAP_2', 'annot' , 'Signature1'))
  
  ggplot(data = tmp, aes(x=annot, y=Signature1))+ geom_boxplot(aes(fill=annot))+
    scale_fill_manual(values=c( "#d62728", "lightpink", "deeppink4"), name=NULL) + 
    labs(x = "", y = "", title = title, fill=title)  +
    theme( plot.title = element_text(size = 10))
}

go_GSEA_df <- as.data.frame(go_GSEA)
oxidative_phosphorylation <- unlist(strsplit(gsub("/", ",", go_GSEA_df$core_enrichment[go_GSEA_df$Description=="oxidative phosphorylation"]),","))
oxidative_phosphorylation_list <- list(oxidative_phosphorylation)

aerobic_respiration <- unlist(strsplit(gsub("/", ",", go_GSEA_df$core_enrichment[go_GSEA_df$Description=="aerobic respiration"]),","))
aerobic_respiration_list <- list(aerobic_respiration)
  
ATP_process <- unlist(strsplit(gsub("/", ",", go_GSEA_df$core_enrichment[go_GSEA_df$Description=="ATP biosynthetic process"]),","))
ATP_process_list <- list(ATP_process)

a <- gsvaplot_cluster(CD14, oxidative_phosphorylation_list, "oxidative phospohorylation")
b <- gsvaplot_cluster(CD14, aerobic_respiration_list, "aerobic respiration")
c <- gsvaplot_cluster(CD14, ATP_process_list, "ATP biosynthetic process")

ggarrange(a,b,c, common.legend =T, nrow = 1)


# C) Boxplots Cell exhaustion per group of patients
gsvaplot_group <- function(Celltype, features_list, title) {
  tmp <- FetchData(
    AddModuleScore(CD14, features=features_list, name = 'Signature', assay='RNA'),
    vars=c('UMAP_1', 'UMAP_2', 'annot' , 'Signature1'))
  
  tmp2 <- FetchData(
    AddModuleScore(CD14, features=features_list, name = 'Signature', assay='RNA'),
    vars=c('UMAP_1', 'UMAP_2', 'groups_new' , 'Signature1'))
  
  tmp$groups <- tmp2$groups_new
  levels(tmp$groups ) <- c("HC","SDC","UDC","pre-ACLF")
   ggplot(data = tmp, aes(x=annot, y=Signature1, fill=factor(groups))) + geom_boxplot()+
    scale_fill_manual(values=c("darkgreen", "gold1", "darkorange1", "firebrick4")) + 
    labs(x = "", y = "", title = title, fill="")  +
     theme( plot.title = element_text(size = 10))
}

a <- gsvaplot_group(CD14, oxidative_phosphorylation_list, "oxidative phospohorylation")
b <- gsvaplot_group(CD14, aerobic_respiration_list, "aerobic respiration")
c <- gsvaplot_group(CD14, ATP_process_list, "ATP biosynthetic process")

ggarrange(a,b,c, common.legend =T, nrow = 1)

# D) Boxplots immunity

inflammatory_response <- unlist(strsplit(gsub("/", ",", go_GSEA_df$core_enrichment[go_GSEA_df$Description=="positive regulation of inflammatory response"]),","))
inflammatory_response_list <- list(inflammatory_response)

bacterium <- unlist(strsplit(gsub("/", ",", go_GSEA_df$core_enrichment[go_GSEA_df$Description=="defense response to bacterium"]),","))
bacterium_list <- list(bacterium)
  
immune <- unlist(strsplit(gsub("/", ",", go_GSEA_df$core_enrichment[go_GSEA_df$Description=="cytokine production involved in immune response"]),","))
immune_list <- list(immune)

a <- gsvaplot_cluster(CD14, inflammatory_response, "positive regulation of inflammatory response")
b <- gsvaplot_cluster(CD14, bacterium_list, "defense response to bacterium")
c <- gsvaplot_cluster(CD14, immune_list, "cytokine production involved in immune response")

ggarrange(a,b,c, common.legend =T, nrow = 1)
