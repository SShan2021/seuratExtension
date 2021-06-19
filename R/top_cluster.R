#' @title top_cluster
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

top_cluster <- function(data){

  final <- data %>% group_by(cluster) %>% top_n(n=1, wt = clusteravg_log2FC)

  return(as.data.frame(final))

}
