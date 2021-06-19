#' @title correctlipa_fisher.test
#'
#' @description Test for clusters between GFP+ and GFP- whether there is signficiantly different/more/less
#' expression of Lipa.
#'
#' @param output A data frame from the output of lipa_expressionmarkers
#' @param counts A data frame from the output of originalcount (gives counts by cluster split by GFP)
#' @param alt For the Fisher's Test: "two.sided" or "greater" or "less" (default is two.sided)
#'
#' @return
#'
#' @examples
#'lipa.40pc_0.2_all <- correctlipa_fisher.test(lipa.40pc_0.2, counts_40pc_0.2)
#'
#' @export
#' @importFrom dplyr "%>%"
#'

correctlipa_fisher.test <- function(output, counts, alt = "two.sided"){

  Lipa_prop = left_join(output, counts) %>%
    mutate(
      GFPneg.n_lipa = round(GFPnegative.pct*GFPnegative),
      GFPpos.n_lipa = round(GFPpositive.pct*GFPpositive)) %>%
    dplyr::select(cluster, GFPpos.n_lipa, GFPneg.n_lipa, GFPpositive, GFPnegative,
                  GFPnegative.pct, GFPpositive.pct, p_val, p_val_adj, clusterlabels)

  pos_sums <- colSums(Lipa_prop[,2:5])
  neg_macrophage <- Lipa_prop %>% filter(str_detect(clusterlabels, "macrophages"))
  neg_sums <- colSums(neg_macrophage[,2:5])

  a_1 <- cbind(pos_sums[1], neg_sums[2])
  a_2 <- cbind(pos_sums[3]-pos_sums[1], neg_sums[4]-neg_sums[2])

  a <- rbind(a_1, a_2)
  rownames(a) <- c("Lipa", "no_Lipa")
  colnames(a) <- c("GFPPos", "GFPNeg")
  return(unlist(fisher.test(a, alternative=alt)[1:3]))

}
