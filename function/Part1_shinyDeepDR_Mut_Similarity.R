#Part1: Nearest neighbors functions-mutation

# load("Datasets/CCLE_22Q2/ccle_22Q2_mutation_cells_1771.RData")
# 
# input.data <- ccle.mut.22Q2[,c(1,20:25)]
# select.sample <- colnames(input.data)[4]
# 
# 
# compaired.table<- "Datasets/Part1_Nearest_Neighbors/gdsc_and_ccle_704_overlapped_cells_mut.RDS"
# 
# system.time({
# a <- shinyDeepDR_Mut_fisher_test(input.data, select.sample, compaired.table = compaired.table)
# })
# 
# 
# #TCGA
# load("Datasets/TCGA/tcga_mutation_samples_9059.RData")
# 
# input.data = tcga.mut[,c(1,3:6)]
# system.time({
#   a <- shinyDeepDR_Mut_fisher_test(input.data = input.data, 
#                                    select.sample = tcga.mut.sampleID[4], 
#                                    compaired.table = tcga.mut)
# })
# 


ShinyDeepDR_Mut_Similarity<- function(input.data, select.sample,
                                       compaired.table, method = c("Fisher","Jaccard")){
  
  if(is.na(select.sample)){
    data = input.data[,c(1,2)]  
    select.sample = colnames(data)[2]
  }else{
    data = input.data[,which(colnames(input.data) %in% c("Gene", select.sample))]
  }
  
  
  #mutation table
  if(class(compaired.table) == "character"){
    compaired.table <- readRDS(compaired.table)
    rownames(compaired.table) <- compaired.table$Gene
  }
  

  if (colnames(compaired.table)[2] == "missing_val") {
    compaired.table <- compaired.table[,-2]
  }
  
  #Matke sure the order of genes in input data and comparied data is the same
  if (sum(data$Gene== compaired.table$Gene) != 18281) {
    data <- data[match(data$Gene, compaired.table$Gene),]
  }
  
  #Fisher
  if(tolower(method) == "fisher"){
  #generate contigency table
  df <- data.frame(DepMapID = colnames(compaired.table)[-1], Fisher_pvalue=NA,
                   #adjPval =NA,
                   stringsAsFactors = F)
  
  tb <- matrix(data = NA, nrow = 2,ncol = 2, dimnames = list(c("User_mut","User_No_mut"),c("compaired_mut","compaired_No_mut")))
  
  
  for (i in 2:ncol(compaired.table)) {
    # #Overlapped
    tb[1,1] <- length(intersect(data$Gene[which(data[,2] == 1)], compaired.table$Gene[which(compaired.table[,i] == 1)]))
    #User mutation only
    tb[1,2] <- length(intersect(data$Gene[which(data[,2] == 1)], compaired.table$Gene[which(compaired.table[,i] == 0)]))
    #Compared data only
    tb[2,1] <-  length(intersect(data$Gene[which(data[,2] == 0)], compaired.table$Gene[which(compaired.table[,i] == 1)]))
    #No mutation
    tb[2,2] <-  length(intersect(data$Gene[which(data[,2] == 0)], compaired.table$Gene[which(compaired.table[,i] == 0)]))

    #Overlapped
    # tb[1,1] <- sum(as.numeric(data[,2])*as.numeric(compaired.table[,i]) == 1)
    # #User mutation only
    # tb[1,2] <- length(intersect(data$Gene[which(data[,2] == 1)], compaired.table$Gene[which(compaired.table[,i] == 0)]))
    # #Compared data only
    # tb[2,1] <-  length(intersect(data$Gene[which(data[,2] == 0)], compaired.table$Gene[which(compaired.table[,i] == 1)]))
    # #No mutation
    # tb[2,2] <-  sum((as.numeric(data[,2])+as.numeric(compaired.table[,i]))==0)
    
    
    df$Fisher_pvalue[which(df$DepMapID == colnames(compaired.table)[i])] <- fisher.test(tb, alternative = "greater")$p.value
    
  }
  }
  
  #Jaccard similarity
  if (tolower(method) == "jaccard") {
    #Jaccard function
    ShinyDeepDR_jaccard <- function(a, b) {
      intersection = length(intersect(a, b))
      union = length(a) + length(b) - intersection
      
      return (intersection/union)
    }
    
    #Generate output table
    df <- data.frame(DepMapID = colnames(compaired.table)[-1], Jaccard=NA,
                     stringsAsFactors = F)
    
    for (i in 2:ncol(compaired.table)) {
      df$Jaccard[i-1] <- ShinyDeepDR_jaccard(a = data$Gene[data[,2] == 1],b = compaired.table$Gene[compaired.table[,i]==1])
    }
    
  }
  
  colnames(df)[2] <- "value"
  cat("Fisher: p-value/ Jaccard: Jaccard index \n")
  return(df)
}
