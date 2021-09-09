#' @title top_cluster (All)
#'
#' @description Gives you the top cluster from cluster_averagelog2FC
#'
#' @param data The output of cluster_averagelog2FC.
#'
#' @return A dataframe with the top clusters labeled.
#'
#' @examples top_cluster(dataset)
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
