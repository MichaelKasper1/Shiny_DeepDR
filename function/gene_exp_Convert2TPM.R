#Date: 20221102
#Function: RPKM or FPKM convert to TPM
#make sure inputData is numeric matrix
#20230215: add if loop for 1 sample

ShinyDeepDR_exp_convert2TPM <- function(inputData){
  
  if (colnames(inputData)[1] == "Gene") {
    rownames(inputData) <- inputData$Gene
    #inputData <- inputData[,-1]
    if (sum(ncol(inputData) -1) == 1) {
      total.rpkm <- sum(as.numeric(inputData[,-1]))
    }else{
    total.rpkm <- colSums(inputData[,-1])}
    
  }else{
    
    if (sum(ncol(inputData) -1) == 1) {
      total.rpkm <- sum(as.numeric(inputData[,-1]))
    }else{
      total.rpkm <- colSums(inputData[,-1])}
    
  }
  
  dt.out <- matrix(data = NA, nrow = nrow(inputData), ncol = ncol(inputData)-1)
  rownames(dt.out) <- inputData$Gene
  colnames(dt.out) <- colnames(inputData)[-1]
  
  if(sum(ncol(inputData) -1) == 1){
    for (i in 1:ncol(dt.out)) {
      
      dt.out[,i] <- as.numeric(inputData[,-1])*10^6/total.rpkm
    }
  }else{
  for (i in 1:ncol(dt.out)) {
    
    dt.out[,i] <- as.numeric(inputData[,i+1])*10^6/total.rpkm[i]
  }
    }

    return(dt.out)
    
}


#a <- ShinyDeepDR_exp_convert2TPM(inputData = inputData[,-1])
