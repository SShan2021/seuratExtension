#' @title signficiant_lipafisher
#'
#' @description Tells you which clusters have signficiant lipa expression.
#'
#' @param dataset A data frame which is the output of correctlipa_fisher.test or lipaall_fisher.test
#'
#' @return A data frame with p-value < 0.05
#'
#' @examples
#'
#'signficiant_lipafisher(cluster1_correctlipafisher)
#'
#' @export
#' @importFrom dplyr "%>%"
#'

signficiant_lipafisher <- function(dataset) {
  return(dataset[dataset$p.value < 0.05,])
}
