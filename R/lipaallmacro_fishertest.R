#' @title lipaallmacro_fishertest
#'
#' @description (LipaTg)Tests for Lipa in macrophages populations in GFP+/GFP-.
#'
#' @param output Output of lipa_expressionmarkers
#' @param counts Output of originalcount
#' @param alt "two.sided", "greater", "less"
#'
#' @return A dataframe of the fisher's test output for Lipa in macrophage
#' subpopulations of GFP+/GFP-.
#'
#' @examples lipa.40pc_0.2 <- lipaallmacro_fishertest(plaque.combined_40pc_0.2,
#' counts = counts_40pc_0.2, alt = "two.sided")
#'
#'
#'
#' @export
#' @importFrom dplyr "%>%"
#'

lipaallmacro_fishertest <- function(output, counts, alt = "two.sided"){

  Lipa_prop = left_join(output, counts) %>%
    mutate(
      GFPneg.n_lipa = round(GFPnegative.pct*GFPnegative),
      GFPpos.n_lipa = round(GFPpositive.pct*GFPpositive)) %>%
    dplyr::select(cluster, GFPpos.n_lipa, GFPneg.n_lipa, GFPpositive, GFPnegative,
                  GFPnegative.pct, GFPpositive.pct, p_val, p_val_adj, clusterlabels)

  Lipa_prop <- Lipa_prop[grep("macrophages", Lipa_prop$clusterlabels),]

  sums <- colSums(Lipa_prop[,2:5])

  a_1 <- rbind(sums[1], sums[3]-sums[1])
  a_2 <- rbind(sums[2], sums[4]-sums[2])
  a <- cbind(a_1,a_2)

  rownames(a) <- c("Lipa", "no_Lipa")
  colnames(a) <- c("GFPPos", "GFPNeg")
  return(unlist(fisher.test(a, alternative=alt)[1:3]))

}
