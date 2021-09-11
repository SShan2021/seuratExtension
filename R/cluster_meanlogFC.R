#' @title cluster_meanlogFC
#'
#' @description (All) Tells you the average meanlogFC per cluster for each celltype
#' (weighted or not weighted )
#'
#' @param dataset The output of identify_cluster(findtop_conservedmarkers())
#' @param weight If you want to have weights, then weight = "weighted"
#'
#' @return A dataframe with the average meanlogFC per cluster for each celltype.
#'
#' @examples cluster_meanlogFC(clust_40pc_0.2_identifymarkers)
#'
#'
#'
#' @export
#' @importFrom dplyr "%>%"
#'

cluster_meanlogFC <- function(dataset, weight = "na"){

  clustercells <- c("all.macrophages", "reslike.macrophages", "inflammatory.macrophages", "trem2high.macrophages",
                    "ifnic.macrophages", "monocyte", "modc.macrophages", "cDC1", "matureDC", "t.cell", "cxcr6.t.cell", "cd8.t.cell",
                    "b.cells", "mast.cells", "granulocytes", "nk.cells", "cell.cycle.related", "smc.cells")

  datalist = list()

  if(weight == "weighted"){
  for(i in 1:length(clustercells)){
    datalist[[i]] <- as.data.frame(dataset %>% filter(celltype == clustercells[i]) %>%
                                     group_by(cluster) %>% summarise(clusteravg_log2FC = mean(mean_logFC*(pct.1-pct.2)), type = clustercells[i]))

  }
  }
  else{
    for(i in 1:length(clustercells)){
      datalist[[i]] <- as.data.frame(dataset %>% filter(celltype == clustercells[i]) %>%
                                       group_by(cluster) %>% summarise(clusteravg_log2FC = mean(mean_logFC), type = clustercells[i]))

    }
  }

  merged_dataframe = do.call(rbind, datalist)

  clustermatrix <- as.data.frame(merged_dataframe %>% group_by_at(vars(cluster:type)) %>%
                                   summarize_all(paste, collapse=","))

  clustermatrix$clusteravg_log2FC <- round(clustermatrix$clusteravg_log2FC,4)

  return(clustermatrix)
}
