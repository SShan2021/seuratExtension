#' @title lipa_fisher.test
#'
#' @description
#'
#' @param dataset
#'
#' @return
#'
#' @examples
#'lipa.30pc_0.4_fisher <- lipa_fisher.test(lipa.30pc_0.4, counts_30pc_0.4)
#'
#'
#' @export
#' @importFrom dplyr "%>%"
#'

lipa_fisher.test <- function(output, counts, alt = "two.sided"){
  Lipa_prop = left_join(output, counts) %>%
    mutate(
      GFPneg.n_lipa = round(GFPnegative.pct*GFPnegative),
      GFPpos.n_lipa = round(GFPpositive.pct*GFPpositive)) %>%
    dplyr::select(cluster, GFPneg.n_lipa, GFPpos.n_lipa, GFPnegative,
                  GFPpositive, GFPnegative.pct, GFPpositive.pct, p_val, p_val_adj, clusterlabels)

  datalist = list()

  for(i in 1:length(Lipa_prop$cluster)){

    lipa_neg <-  unlist(Lipa_prop[i,c("GFPpositive", "GFPnegative")]) -
      unlist(Lipa_prop[i,c("GFPpos.n_lipa", "GFPneg.n_lipa")])


    datalist[[i]] <- rbind(unlist(Lipa_prop[i,c("GFPpos.n_lipa", "GFPneg.n_lipa")]),
                           lipa_neg)

    rownames(datalist[[i]]) <- c("Lipa", "no_Lipa")
    colnames(datalist[[i]]) <- c("GFPPos", "GFPNeg")
  }

  newlist = list()
  for(i in 1:length(datalist)){
    newlist[[i]] <- fisher.test(datalist[[i]], alternative=alt)
  }

  merged <- as.data.frame(do.call(rbind, newlist))
  merged <- merged[,1:3]
  merged$cluster <- Lipa_prop$cluster
  merged$clusterlabels <- Lipa_prop$clusterlabels

  return(merged)

}
