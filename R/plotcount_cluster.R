#' @title plotcount_cluster
#'
#' @description
#'
#' @param dataset
#'
#' @return
#'
#' @examples
#'
#'
#'
#' @export
#' @importFrom dplyr "%>%"
#'

plotcount_cluster <- function (dataset) {
  dataset@meta.data %>% dplyr::count(seurat_clusters, clusterlabels, GFP) %>%
    ggplot(aes(x = seurat_clusters, y = n, fill = GFP)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label=n), position=position_dodge(width=0.9), size=3) +
    labs(title = "Distribution of GFP+/GFP- Cells", x = "Cluster ID")
}
