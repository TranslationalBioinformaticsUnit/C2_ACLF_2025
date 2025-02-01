library(seurat)
library(ggplot)
library(dplyr)



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

color_palette <-  (c( "#d62728", "lightpink", "deeppink4"))

proportions_df_cells <- find_proportions_df(
  CD14,
  x = "groups",
  fill = "annotation"
)

stacked_barplot_cells <- plot_stacked_barplot(
  proportions_df_cells,
  x = "groups",
  fill = "annotation",
  colors = color_palette
)
