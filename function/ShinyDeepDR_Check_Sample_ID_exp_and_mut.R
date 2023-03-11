#Function: check sample name of gene mutation data & gene expression data

# setwd("~/OneDrive - University of Pittsburgh/UPitt_project/shinyDeepDR_work/")
# 
# #example
# exp <- read.delim("Datasets/example_data_for_test_function/gene_exp/gene_exp_fpkm_5samples_mix.txt", sep = "\t", stringsAsFactors = F)
# exp.data.values = "FPKM"
# mut <- read.delim("Datasets/example_data_for_test_function/gene_mut/gene_mut_matrix_5samples_mix.txt",sep = "\t", stringsAsFactors = F)
# mut.data.format = "MAF"
# a <- ShinyDeepDR_CheckSampleNames(exp,mut ,mut.data.format = "TXT")

ShinyDeepDR_CheckSampleNames <- function(exp, mut, 
                                         mut.data.format = c("TXT", "MAF","LIST")){
  
  #Gene expression
  colnames(exp)[1] <- "Gene"
  sampleID.exp <- colnames(exp)[-1]
  
  #Mutation
  #TXT format
  if (toupper(mut.data.format) == "TXT") {
    colnames(mut)[1] <- "Gene"
    sampleID.mut <- colnames(mut)[-1]
  }
  
  #MAF format
  if (toupper(mut.data.format) == "MAF") {
    if (sum(colnames(mut) %in% c("Hugo_Symbol" ,"Variant_Classification","Tumor_Sample_Barcode"  )) == 3) {
      sampleID.mut <- unique(mut$Tumor_Sample_Barcode)
      
    }else{
      idx <- which( c("Hugo_Symbol" ,"Variant_Classification","Tumor_Sample_Barcode"  ) %in% colnames(mut))
 
      cat( paste(c("Hugo_Symbol" ,"Variant_Classification","Tumor_Sample_Barcode")[-idx],collapse = " & "),"column missing!  \n")
    }
    
  }
  
  #List
  if (toupper(mut.data.format) == "LIST") {
    
    if(length(sampleID.exp) == 1){
    sampleID.mut <- sampleID.exp
    
  }else{
    cat("The mutated gene list is for a single sample only!")
  }
  }
  
  
  #Check both sampleID
      commonID <-intersect(sampleID.exp,sampleID.mut)
      
      #The sample names on both files do not match
      if (length(commonID) == 0) {
        cat("The sample names on both files do not match!! \n")
      }
      
      ##Part of sample names om both files do not match
      if (length(commonID) != length(sampleID.exp) | length(commonID) != length(sampleID.mut) ) {
        cat("Part of sample names on both files do not match!! \n")
      }

return(commonID)
  
}
