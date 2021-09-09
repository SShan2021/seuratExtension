#' @title significant_lipafisher
#'
#' @description (LipaTg) Gives you which cluster have significant
#'
#' @param dataset The output of correctlipa_fisher.test
#'
#' @return A dataframe with the significant Lipa clusters
#'
#' @examples lipa.40pc_0.2_sig <- significant_lipafisher(lipa.40pc_0.2_fisher)
#'

#' @export
#' @importFrom dplyr "%>%"
#'

significant_lipafisher <- function(dataset) {
  return(dataset[dataset$p.value < 0.05,])
}

