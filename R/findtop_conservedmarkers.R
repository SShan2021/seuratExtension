#' @title findtop_conservedmarkers
#'
#' @description (LipaTg)Finds the top 20 conserved markers to be put into
#' identify_cluster.
#'
#' @param datalist A data frame of conserved markers
#'
#' @return A dataframe of the top 20 conserved markers
#'
#' @examples clust_40pc_0.4_identifymarkers <- identify_cluster(findtop_conservedmarkers(clust_40pc_0.4_conservedmarkers))
#'
#'
#' @export
#' @importFrom dplyr "%>%"
#'

findtop_conservedmarkers <- function(datalist){

  for(i in 1:length(datalist)){
    if(!is.null(datalist[[i]])){
      datalist[[i]]$gene <- rownames(datalist[[i]])
      datalist[[i]]$cluster <- i-1
    }
  }

  datalist[sapply(datalist, is.null)] <- NULL
  newlist = list()

  for(i in 1:length(datalist)){
    if(!is.null(datalist[[i]]$GFPpositive_avg_log2FC) && !is.null(datalist[[i]]$GFPnegative_avg_log2FC)){
      newlist[[i]] <- as.data.frame(datalist[[i]]) %>%
        dplyr::mutate(mean_logFC = log((exp(GFPpositive_avg_log2FC) + exp(GFPnegative_avg_log2FC))/2)) %>%
        arrange(desc(mean_logFC)) %>% dplyr::select(gene, cluster, mean_logFC) %>% top_n(n=20, wt = mean_logFC)
    }
  }
  merged_dataframe = do.call(rbind, newlist)

  return(merged_dataframe)
}
