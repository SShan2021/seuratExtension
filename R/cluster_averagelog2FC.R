#' @title cluster_averagelog2FC
#'
#' @description (All)Goes through given top markers and outputs the avg_log2FC per cluster for each cell type
#'
#' @param dataset The data frame output of identify_cluster.R
#'
#' @return A data frame object that has avg_log2FC per cluster for each
#'
#' @examples
#' cluster2_avglog2FC <- cluster_averagelog2FC(dataset = output_table, weight = "weighted")
#'
#'
#' @export
#' @importFrom dplyr "%>%"
#'

cluster_averagelog2FC <- function(dataset, weight = "na"){

  clustercells <- c("all.macrophages", "reslike.macrophages", "inflammatory.macrophages", "trem2high.macrophages",
                    "ifnic.macrophages", "monocyte", "modc.macrophages", "cDC1", "matureDC", "t.cell", "cxcr6.t.cell", "cd8.t.cell",
                    "b.cells", "mast.cells", "granulocytes", "nk.cells", "cell.cycle.related", "smc.cells")

  datalist = list()
  if(weight == "na"){
  for(i in 1:length(clustercells)){
    datalist[[i]] <- as.data.frame(dataset %>% filter(celltype == clustercells[i]) %>%
                                     group_by(cluster) %>% summarise(clusteravg_log2FC = mean(avg_log2FC), type = clustercells[i]))

  }}
  else{
    for(i in 1:length(clustercells)){
      datalist[[i]] <- as.data.frame(dataset %>% filter(celltype == clustercells[i]) %>%
                                       group_by(cluster) %>% summarise(clusteravg_log2FC = mean(avg_log2FC*(pct.1-pct.2)), type = clustercells[i]))

    }
  }

  merged_dataframe = do.call(rbind, datalist)

  clustermatrix <- as.data.frame(merged_dataframe %>% group_by_at(vars(cluster:type)) %>%
                                   summarize_all(paste, collapse=","))

  clustermatrix$clusteravg_log2FC <- round(clustermatrix$clusteravg_log2FC,4)

  return(clustermatrix)
}
