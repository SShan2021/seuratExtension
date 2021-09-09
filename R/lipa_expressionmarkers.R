#' @title lipa_expressionmarkers (LipaTg)
#'
#' @description Returns the expression per cluster of the feature you
#' want expressed sorted by avg2logFC
#'
#' @param seuratobject A seurat object.
#' @param clusternumber Output of top_cluster.
#' @param feature Name of the feature you want to find
#'
#' @return A dataframe of the avg2logFC of the feature you want expressed
#' in each cluster.
#'
#' @examples lipa.40pc_0.2 <- lipa_expressionmarkers(plaque.combined_40pc_0.2,clust_40pc_0.2,
#' clusternumber = clust_40pc_0.2, feature = "Lipa")
#'
#'
#'
#' @export
#' @importFrom dplyr "%>%"
#'

lipa_expressionmarkers <- function(seuratobject, clusternumber, feature){
  ## Assign identify based on GFP+/GFP- group and cluster membership
  seuratobject$celltype <- paste(Idents(seuratobject), seuratobject$GFP, sep = "_")
  Idents(seuratobject) <- "celltype"

  datalist = list()

  for (i in 0:length(clusternumber[,1])) {

    skip_to_next <- FALSE

    tryCatch(
      datalist[[i]] <- FindMarkers(seuratobject, paste0(i-1, "_GFPpositive"),
                                   ident.2 = paste0(i-1, "_GFPnegative"), features = feature,
                                   logfc.threshold = 0, min.diff.pct = 0, min.pct = 0,
                                   min.cells.feature = 0, min.cells.group = 0),
      error = function(e)
      { skip_to_next <<- TRUE
      message('Caught an error!')
      print(e)})

    if(skip_to_next) { next }
  }

  for(i in 1:length(datalist)){
    if(!is.null(datalist[[i]])){
      datalist[[i]]$gene <- rownames(datalist[[i]])
      datalist[[i]]$cluster <- i-1
      datalist[[i]]$clusterlabels <- clusternumber$type[i]
    }
  }

  datalist[sapply(datalist, is.null)] <- NULL

  merged <- do.call(rbind, datalist)
  rownames(merged) <- 1:length(merged$cluster)

  merged$cluster <- as.factor(merged$cluster)

  return(merged %>%
           dplyr::select(gene, cluster, clusterlabels, avg_log2FC, GFPpositive.pct = pct.1, GFPnegative.pct = pct.2, p_val, p_val_adj))

}
