#' @title celltypecount
#'
#' @description (LipaTg)Outputs the counts in the cluster by celltype given by
#' identify_cluster
#'
#' @param dataset A seurat object.
#' @param label Original datatype (ie. GFP+/GFP-)
#'
#' @return A dataframe of the counts in the cluster by celltype.
#'
#' @examples cellcounts_40pc_0.2 <- celltypecount(plaque.combined_40pc_0.2)
#'
#'
#'
#' @export
#' @importFrom dplyr "%>%"
#'

celltypecount <- function(dataset, label){
  counts <- dataset@meta.data %>% dplyr::count(clusterlabels, GFP) %>% spread(clusterlabels, n)
  return(counts)
}
