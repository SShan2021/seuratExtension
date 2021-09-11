#' @title plotcount_cluster
#'
#' @description (LipaTg)Plots the distribution of cells by cluster, split by
#' original datatype.
#'
#' @param dataset A seurat object.
#' @param label Original datatype (ie. GFP+/GFP-)
#'
#' @return A barplot of the counts in each cluster by celltype.
#'
#' @examples plotcounts_40pc_0.2 <- plotcount_cluster(plaque.combined_40pc_0.2)
#'
#'
#'
#' @export
#' @importFrom dplyr "%>%"
#'

plotcount_cluster <- function (dataset, GFP) {

  a <- dataset@meta.data %>% dplyr::count(seurat_clusters, clusterlabels, GFP) %>%
    group_split(GFP)

  b <- a[[1]] %>% mutate(sum = sum(n)) %>% mutate (prop = n/sum)
  c <- a[[2]] %>% mutate(sum = sum(n)) %>% mutate (prop = n/sum)

  rbind(b,c) %>%
    ggplot(aes(x = seurat_clusters, y = prop, fill = GFP)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label=round(prop,2)), position=position_dodge(width=0.9), size=3) +
    labs(title = "Distribution of GFP+/- cells", x = "Cluster ID") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
