#' @title originalcount
#'
#' @description (All)Outputs the counts in the cluster by original datatype.
#'
#' @param dataset A seurat object.
#'
#' @return A dataframe of the counts in the cluster by original datatype
#'
#' @examples counts_40pc_0.2 <- originalcount(plaque.combined_40pc_0.2)
#' would return counts split by GFP+/GFP-
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
