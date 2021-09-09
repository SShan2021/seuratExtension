#' @title findall_conservedmarkers (All)
#'
#' @description Finds all the conserved markers for dataset automatically skipping
#' datasets which are too small
#'
#' @param seuratobject The data frame output of identify_cluster.R
#' @param clusternumber The number of clusters
#' @param grouping Grouping Variable (ie. GFP)
#'
#' @return A data frame object that has avg_log2FC per cluster for each
#'
#' @examples
#' findall_conservedmarkers(harmony_30pc_20param_0.2, clusternumber = 9, grouping = "GFP")
#'
#'
#' @export
#' @importFrom dplyr "%>%"
#'

findall_conservedmarkers <- function(seuratobject, clusternumber, grouping = "sample"){
  datalist = list()

  for (i in 1:clusternumber) {

    skip_to_next <- FALSE

    tryCatch(
      datalist[[i]] <- FindConservedMarkers(seuratobject,
                                            ident.1 = i-1, grouping.var = grouping),
      error = function(e)
      { skip_to_next <<- TRUE
      message('Caught an error!')
      print(e)})

    if(skip_to_next) { next }
  }

  return(datalist)

}
