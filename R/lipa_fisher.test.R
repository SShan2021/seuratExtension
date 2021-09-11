#' @title lipa_fisher.test
#'
#' @description (LipaTg) Fisher’s test on each of the clusters 
#'
#' @param output Output of lipa_expressionmarkers
#' @param counts Output of originalcount
#' @param alt "two.sided", "greater", "less"
#'
#' @return A dataframe of the fisher's test output for Lipa in each
#' subpopulations of GFP+/GFP-.
#'
#' @examples lipa_fisher.test(lipa.40pc_0.2, counts_40pc_0.2)
#'
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
    dplyr::select(cluster, avg_log2FC, GFPneg.n_lipa, GFPpos.n_lipa, GFPnegative,
                  GFPpositive, GFPnegative.pct, GFPpositive.pct, p_val, p_val_adj, clusterlabels)
  
  Lipa_prop <- na.omit(Lipa_prop)
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
  merged$avg_log2FC <- Lipa_prop$avg_log2FC
  
  return(merged)
  
}