#' @title rename_indents
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
