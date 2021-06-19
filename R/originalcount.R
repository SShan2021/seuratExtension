#' @title originalcount
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

originalcount <- function(dataset){
  counts <- dataset@meta.data %>% dplyr::count(seurat_clusters, clusterlabels, GFP) %>% spread(GFP, n)
  counts <- counts %>% dplyr::rename(cluster = seurat_clusters)
  return(counts)
}
