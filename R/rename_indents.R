#' @title rename_indents
#'
#' @description (All)Relabels the clusters automatically for graphing.
#'
#' @param seuratobject The seurat object and the outuput of top_cluster.
#' @param labels The output of top_cluster
#'
#' @return N/a
#'
#' @examples plaque.combined_40pc_0.6 <- rename_indents(plaque.combined_40pc_0.6,
top_cluster(plaque.combined_40pc_0.6clusteraverages))
#'
#'
#'
#' @export
#' @importFrom dplyr "%>%"
#'

rename_indents <- function(seuratobject, labels){


  data <- data.frame(seurat_clusters = seuratobject@meta.data$seurat_clusters)
  clusternumber <- as.data.frame(table(data))
  colnames(clusternumber) <- c("cluster", "freq")

  together <- merge(clusternumber, labels, all = TRUE) %>% arrange(cluster)

  values <- together$type
  data$clusterlabels <- values[data$seurat_clusters]
  seuratobject@meta.data$clusterlabels <- data$clusterlabels

  return(seuratobject)
}
