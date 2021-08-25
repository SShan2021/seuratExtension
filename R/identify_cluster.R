#' @title identify_cluster
#'
#' @description Goes through given top markers and extract those that match the
#' cluster identifiers.
#'
#' @param topmarkers A dataset of topmarkers
#'
#' @return A data frame object that has all the markers and corresponding
#' clusters and the cluster type for each gene.
#'
#' @examples
#' output_table <- identify_cluster(topmarkers, upper = "upper", top20 = "top20", writecsv = "writecsv")
#'
#'
#' @export
#' @importFrom dplyr "%>%"
#'


identify_cluster <- function(topmarkers, upper = "na", top20 = "na", writecsv = "na"){

  source('http://www.dr-spiess.de/scripts/cbind.na.R')

  all.macrophages <- c("Adgre1", "Csf1r", "Fcgr1", "CD64", "Cd68", "Cd14", "Mafb", "Ly6c2")
  reslike.macrophages <- c("Csf1r", "Pf4", "Txnip", "F13a1", "Sepp1", "Lyve1", "Timp2", "Lyz2", "Folr2", "Ccl9", "Gas6", "CD206", "Mrc1")
  inflammatory.macrophages <- c("Cxcl2", "Cd14", "Cebpb", "Nlrp3", "Il1b", "Nfkbiz", "Egr1", "Zpf36", "Ccl2", "Tnf", "Ier3", "Ccl3", "Tlr2", "Nfkbid")
  trem2high.macrophages <- c("Trem2", "Cd9", "Ctsd", "Spp1", "Lgals3", "Atox1", "Ctsb", "Ctsz")
  ifnic.macrophages <- c("Isg15", "Irf7")
  monocyte <- c("F10", "Lilrb4a", "Ly6c2", "Ccr2", "Csf1r")
  modc.macrophages <- c("Cd209a", "Flt3", "Itgb7", "Napsa", "Ifi30", "Syngr2", "Cd74", "H2-Eb1")
  cDC1 <- c("Wdfy4", "Cd24", "Cd74")
  matureDC <- c("Ccr7", "Fscn1")
  t.cell <- c("Lck", "Cd3d", "Rag1", "Cd8a")
  cxcr6.t.cell <- c("Cxcr6", "Il7r", "Icos", "Cd3g")
  cd8.t.cell <- c("Cd8b1", "Cd3g", "Nkg7", "Cd8a")
  b.cells <- c("Cd79a", "Mzb1", "Cd79b", "Ly6d", "Ebf1", "Plac8")
  mast.cells <- c("Furin", "Il1rl1", "Calca")
  granulocytes <- c("Ngp", "Camp", "S100A8", "S100A9")
  nk.cells <- c("Klrb1c", "Gzma", "Klrc1", "Ncr1")
  cell.cycle.related <- c("Ccna2", "Cdk1", "Cdk4")
  smc.cells <- c("Acta2", "Cnn1", "Tagln", "Myh11")
  endothelial.cells <- c("Cd31", "Cd34", "Cd45", "Cd54", "Lyve1",
                         "Tek", "Cd106")


  if(upper == "upper"){
    all.macrophages <- toupper(all.macrophages)
    reslike.macrophages <- toupper(reslike.macrophages)
    inflammatory.macrophages <- toupper(inflammatory.macrophages)
    trem2high.macrophages <- toupper(trem2high.macrophages)
    ifnic.macrophages <- toupper(ifnic.macrophages)
    monocyte <- toupper(monocyte)
    modc.macrophages <- toupper(modc.macrophages)
    cDC1 <- toupper(cDC1)
    matureDC <- toupper(matureDC)
    t.cell <- toupper(t.cell)
    cxcr6.t.cell <- toupper(cxcr6.t.cell)
    cd8.t.cell <- toupper(cd8.t.cell)
    b.cells <- toupper(b.cells)
    mast.cells <- toupper(mast.cells)
    granulocytes <- toupper(granulocytes)
    nk.cells <- toupper(nk.cells)
    cell.cycle.related <- toupper(cell.cycle.related)
    endothelial.cells <- toupper(endothelial.cells)
  }

  clustermarkers_list <- c(all.macrophages, reslike.macrophages, inflammatory.macrophages, trem2high.macrophages,
                           ifnic.macrophages, monocyte, modc.macrophages, cDC1, matureDC, t.cell, cxcr6.t.cell, cd8.t.cell,
                           b.cells, mast.cells, granulocytes, nk.cells, cell.cycle.related, smc.cells,
                           endothelial.cells)

  clustermarkers <- as.data.frame(cbind.na(all.macrophages, reslike.macrophages, inflammatory.macrophages, trem2high.macrophages,
                                           ifnic.macrophages, monocyte, modc.macrophages, cDC1, matureDC, t.cell, cxcr6.t.cell,
                                           cd8.t.cell,b.cells, mast.cells, granulocytes, nk.cells, cell.cycle.related, smc.cells,
                                           endothelial.cells))

  if(top20 == "top20"){
    topmarkers <- as.data.frame(topmarkers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))
  }

  datalist = list()
  for(i in 1:length(clustermarkers_list)){
    datalist[[i]] <- topmarkers[topmarkers$gene == clustermarkers_list[i],]
  }

  merged_dataframe = do.call(rbind, datalist)

  merged_dataframe$celltype <- "other"

  for(i in 1:length(colnames(clustermarkers))){
    for(j in 1:length(na.omit(clustermarkers[,i]))){
      for(k in 1:length(merged_dataframe$gene)){
        if(merged_dataframe$gene[k] == clustermarkers[j,i]){
          merged_dataframe$celltype[k] <- colnames(clustermarkers)[i]
        }
      }
    }

  }

  rownames(merged_dataframe) <- 1:length(merged_dataframe$gene)

  if(!is.null(merged_dataframe$avg_log2FC)){
    merged_dataframe <- merged_dataframe %>% arrange(desc(avg_log2FC)) %>% group_by(cluster) %>%
      summarise(p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, gene, celltype) %>%
      dplyr::select(p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, gene, celltype)

    write.csv(topmarkers %>% group_by(cluster) %>% top_n(n=20, wt = avg_log2FC), "topDEgenes.csv")
  }
  else{
    merged_dataframe <- merged_dataframe %>% arrange(desc(mean_logFC)) %>% group_by(cluster) %>%
      summarise(gene, mean_logFC, celltype) %>%
      dplyr::select(mean_logFC, cluster, gene, celltype)

  }

  if(writecsv == "writecsv"){
    write.csv(topmarkers, "topDEconservedgenes.csv")
  }

  merged_dataframe <- as.data.frame(merged_dataframe)

  merged_dataframe <- unique(merged_dataframe)

  return(merged_dataframe)

}
