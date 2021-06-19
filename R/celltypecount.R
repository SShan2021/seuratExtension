#' @title celltypecount
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

celltypecount <- function(dataset){
  counts <- dataset@meta.data %>% dplyr::count(clusterlabels, GFP) %>% spread(clusterlabels, n)
  return(counts)
}
