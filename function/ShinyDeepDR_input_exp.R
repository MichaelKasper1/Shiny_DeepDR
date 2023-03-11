#Date: 20221102
#Function: check gene expression data
#20230215 
#a. add command line to check order of gene symbols
#b. Filled missing values

#example:
#gene expression
# inputData <-  read.delim("example_data_for_test_function/gene_exp/gene_exp_fpkm_txt", sep = "\t", stringsAsFactors = F)

ShinyDeepDR_input_exp <- function(inputData,data.values = c("TPM", "RPKM","FPKM"),selectID=NA,
                                  log2_transform = c(1, 2),
                                  alias.table = "../ccle_exp_and_mut_with_gene_alias.RData"){
  library(methods)
  load("data/ccle_exp_and_mut_index_missingVal.RData")
  load("data/CCLE_22Q2/ccle_22Q2_expression_cells_1406.RData")
  source("function/CheckGeneSymbol.R")
  source("function/gene_exp_Convert2TPM.R")
  
#created S4 object
  setClass(Class = "input_exp",
           representation = representation(inputData = "data.frame",
                                           inputData.type = "character",
                                           sample_ID = "character",
                                           genelist = "character",
                                           log2_transform = "character",
                                           genelist.update = "data.frame",
                                           inputData.rpkm2tpm = "data.frame",
                                           exp_final = "data.frame"))
  
  colnames(inputData)[1] <- "Gene"
  inputData$Gene <- toupper(inputData$Gene)
  
  if (!is.na(selectID)) {
    inputData <- inputData[,c("Gene", selectID)]
  }
  
  #TPM
  if (data.values == "TPM") {
  if (log2_transform == 1) {
    log2_transform = "YES"
  }else{
    Gene = inputData$Gene
    sampleName <- colnames(inputData)[-1]
    inputData =as.data.frame( cbind(Gene, log2(inputData[,-1]+1)))
    colnames(inputData)[-1] <- sampleName
    log2_transform = "No"
  }
  
  exp.out <- new("input_exp",
                   inputData = inputData,
                   inputData.type = data.values,
                   sample_ID = colnames(inputData)[-1],
                   log2_transform = log2_transform ,
                   genelist.update = ShinyDeepDR_CheckGeneSymbol(inputData = inputData,data.type = "exp",alias.table = alias.table),
                   genelist = inputData$Gene
                   )
  
  #check genes of input data
  if(sum(toupper(inputData$Gene) %in% ccle.exp.index.missingVal$Gene) != nrow(ccle.exp.index.missingVal)){
    idx <- which(!is.na(match(ccle.exp.index.missingVal$Gene,inputData$Gene)))
    tb1 <- ccle.exp.index.missingVal
    tb1[idx,2] <- inputData[,2] 
    colnames(tb1)[2] <- colnames(inputData)[2]
    exp.out@exp_final <- tb1
    
  }else{
    exp.out@exp_final <- exp.out@inputData 
  }
  
  exp.out@exp_final$Gene <-exp.out@genelist.update$gene.symbol.update
  exp.out@exp_final <- exp.out@inputData[match(ccle.exp.22Q2$Gene,exp.out@inputData$Gene),]
  identical(ccle.exp.index.missingVal$Gene, exp.out@exp_final$Gene)
  
  #Filling missing values
  if(sum(is.na(exp.out@exp_final[,2]))>0){
    idx <- which(is.na(exp.out@exp_final[,2]))
    exp.out@exp_final[idx,2] <- ccle.exp.index.missingVal$missing_val[idx]
  }
  
  
}
  
  #convert rpkm or fpkm to tpm
  if (data.values != "TPM") {
    if (log2_transform == 1) {
      log2_transform = "YES"
      Gene = inputData$Gene
      sampleName <- colnames(inputData)[-1]
      inputData =as.data.frame( cbind(Gene, (2^inputData[,-1])-1))
      colnames(inputData)[-1] <- sampleName
      inputData[,-1] <- as.numeric(inputData[,-1] )
      
    }else{
      log2_transform = "No"
    }
    
    exp.out <- new("input_exp",
                   inputData = inputData,
                   inputData.type = data.values,
                   sample_ID = colnames(inputData)[-1],
                   log2_transform = log2_transform ,
                   genelist.update = ShinyDeepDR_CheckGeneSymbol(inputData = inputData,data.type = "exp",alias.table = alias.table),
                   genelist = inputData$Gene
    )
    
    
    rpkm2tpm <- ShinyDeepDR_exp_convert2TPM(inputData = exp.out@inputData)
    colnames(rpkm2tpm ) <- exp.out@sample_ID
    rownames(rpkm2tpm ) <- exp.out@genelist.update$gene.symbol.update
    exp.out@inputData.rpkm2tpm  <-as.data.frame( log2(rpkm2tpm +1) )
    exp.out@exp_final <-  exp.out@inputData.rpkm2tpm  %>% add_column(Gene =  exp.out@genelist.update$gene.symbol.update,
                                                                     .before = colnames(exp.out@inputData.rpkm2tpm )[1] )
    
    #check genes of input data
    if(sum(toupper(exp.out@exp_final$Gene) %in% ccle.exp.index.missingVal$Gene) != nrow(ccle.exp.index.missingVal)){
      idx <- which(!is.na(match(ccle.exp.index.missingVal$Gene,exp.out@exp_final$Gene)))
      tb1 <- ccle.exp.index.missingVal
      tb1[idx,2] <- exp.out@exp_final[,2] 
      colnames(tb1)[2] <- colnames(inputData)[2]
      exp.out@exp_final <- tb1
      
    }
    
    exp.out@exp_final <- exp.out@exp_final[match(ccle.exp.index.missingVal$Gene,exp.out@exp_final$Gene),]
    identical(ccle.exp.index.missingVal$Gene, exp.out@exp_final$Gene)
    
    #Filling missing values
    if(sum(is.na(exp.out@exp_final[,2]))>0){
      idx <- which(is.na(exp.out@exp_final[,2]))
      exp.out@exp_final[idx,2] <- ccle.exp.index.missingVal$missing_val[idx]
    }
    
  }
  
  rownames(exp.out@exp_final) <- exp.out@exp_final$Gene
  return(exp.out)
}

