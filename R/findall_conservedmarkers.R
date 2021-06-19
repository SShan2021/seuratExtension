#' @title findall_conservedmarkers
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

findall_conservedmarkers <- function(seuratobject, clusternumber){
  datalist = list()

  for (i in 1:length(clusternumber)) {

    skip_to_next <- FALSE

    tryCatch(
      datalist[[i]] <- FindConservedMarkers(seuratobject,
                                            ident.1 = clusternumber[i], grouping.var = "GFP"),
      error = function(e)
      { skip_to_next <<- TRUE
      message('Caught an error!')
      print(e)})

    if(skip_to_next) { next }
  }

  return(datalist)

}
