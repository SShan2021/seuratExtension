#' @title cluster_meanlogFC
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

cluster_meanlogFC <- function(dataset){
  #  t1 <- as.data.frame(dataset %>% filter(celltype == "all.macrophages") %>%
  #    group_by(cluster) %>% summarise(mean(avg_log2FC)))


  clustercells <- c("all.macrophages", "reslike.macrophages", "inflammatory.macrophages", "trem2high.macrophages",
                    "ifnic.macrophages", "monocyte", "modc.macrophages", "cDC1", "matureDC", "t.cell", "cxcr6.t.cell", "cd8.t.cell",
                    "b.cells", "mast.cells", "granulocytes", "nk.cells", "cell.cycle.related", "smc.cells")

  datalist = list()
  for(i in 1:length(clustercells)){
    datalist[[i]] <- as.data.frame(dataset %>% filter(celltype == clustercells[i]) %>%
                                     group_by(cluster) %>% summarise(clusteravg_log2FC = mean(mean_logFC), type = clustercells[i]))

  }

  merged_dataframe = do.call(rbind, datalist)

  clustermatrix <- as.data.frame(merged_dataframe %>% group_by_at(vars(cluster:type)) %>%
                                   summarize_all(paste, collapse=","))

  clustermatrix$clusteravg_log2FC <- round(clustermatrix$clusteravg_log2FC,4)

  return(clustermatrix)
}
