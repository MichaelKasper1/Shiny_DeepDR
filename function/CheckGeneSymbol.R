#'Check DepOIs in the fingerprint database


#inputData <- read.delim("example_data_for_test_function/gene_exp/gene_exp_fpkm_txt", sep = "\t", stringsAsFactors = F)

ShinyDeepDR_CheckGeneSymbol <- function(inputData, data.type = c("exp", "mut"),
                            alias.table = "ccle_exp_and_mut_with_gene_alias.RData"){

  library(stringr)
  library(dplyr)
  library(tibble)


  load(alias.table)
  colnames(inputData)[1] <- "Gene"
  geneLists <- toupper(inputData$Gene)

  if(data.type == "exp"){
    check.table <- as.data.frame(apply(ccle.exp.ensembl_gene_alias,MARGIN = 2, FUN = function(x) return(toupper(x))))
    #rownames(check.table) <- check.table$Gene
  }

  if(data.type == "mut"){
    check.table <- as.data.frame(apply(ccle.mut.ensembl_gene_alias,MARGIN = 2, FUN = function(x) return(toupper(x))))
    #rownames(check.table) <- check.table$Gene
  }


  ##identify gene id format
  if (sum(which(substr(geneLists,1,4) == "ENSG")) == 0) {
    #gene symbol only
    gene.symbol <- geneLists[which(substr(geneLists,1,4) != "ENSG")]
    out.table <- data.frame(gene.symbol, gene.symbol.update = NA, stringsAsFactors = F)

    idx <- match(out.table$gene.symbol, check.table$Gene)
    out.table$gene.symbol.update[!is.na(idx)] <-  check.table$Gene[idx[!is.na(idx)]]

    for (i in which(is.na(out.table$gene.symbol.update))) {
      tb1 <- which(check.table == gene.symbol[i], arr.ind = T)
      if(nrow(tb1) == 1){
        out.table$gene.symbol.update[i] <- check.table$Gene[tb1[1]]
      }

    }

  } else if(sum(substr(geneLists,1,4) == "ENSG") == length(geneLists)){

    #ensembl id only
    ensembl <- geneLists[which(substr(geneLists,1,4) == "ENSG")]
    out.table <- data.frame(ensembl, ensembl.update = NA, stringsAsFactors = F)

    idx <- match(out.table$ensembl, check.table$Gene)
    out.table$ensembl.update[!is.na(idx)] <-  check.table$Gene[idx[!is.na(idx)]]

    for (i in which(is.na(out.table$ensembl.update))) {
      tb1 <- which(check.table == ensembl[i], arr.ind = T)
      if(nrow(tb1) == 1){
        out.table$ensembl.update[i] <- check.table$Gene[tb1[1]]
      }
    }


  }else{

  #mix
    gene.symbol <- geneLists[which(substr(igeneLists,1,4) != "ENSG")]
    ensembl <- geneLists[which(substr(geneLists,1,4) == "ENSG")]

    gene.symbol.out <- data.frame(gene.symbol, gene.symbol.update = NA, stringsAsFactors = F)

    idx <- match(out.table$gene.symbol, check.table$Gene)
    out.table$gene.symbol.update[!is.na(idx)] <-  check.table$Gene[idx[!is.na(idx)]]

    for (i in which(is.na(out.table$gene.symbol.update))) {
      tb1 <- which(check.table == gene.symbol[i], arr.ind = T)
      if(nrow(tb1) == 1){
        gene.symbol.out$gene.symbol.update[i] <- check.table$Gene[tb1[1]]
      }

    }

    ensembl.out <- data.frame(ensembl, ensembl.update = NA, stringsAsFactors = F)

    idx <- match(out.table$ensembl, check.table$Gene)
    out.table$ensembl.update[!is.na(idx)] <-  check.table$Gene[idx[!is.na(idx)]]

    for (i in which(is.na(out.table$ensembl.update))) {
      tb1 <- which(check.table == ensembl[i], arr.ind = T)
      if(nrow(tb1) == 1){
        ensembl.out$ensembl.update[i] <- check.table$Gene[tb1[1]]
      }
    }

    out.table <- rbind(gene.symbol.out,ensembl.out)
  }


  return(out.table)
}


