#Date: 20221102
#Function: check gene mutation data

#20230215
#a. remove variant types: "Splice_Site"

# setwd("~/OneDrive - University of Pittsburgh/UPitt_project/shinyDeepDR_work/")
# 
# #example:
# #gene mutation
# #mutation matrix
# inputData <- read.delim("Datasets/example_data_for_test_function/gene_mut/gene_mut_matrix_5samples.txt", sep = "\t", stringsAsFactors = F)
# a <- ShinyDeepDR_input_mut(inputData = inputData,data.format = "TXT")
# 
# #MAF -->DepMap_ID renames to Tumor_Sample_Barcode
# inputData <-  read.delim("Datasets/example_data_for_test_function/gene_mut/data_mutations_mskcc_maf.txt", sep = "\t", stringsAsFactors = F)
# a <- ShinyDeepDR_input_mut(inputData = inputData,data.format = "MAF",alias.table = "ccle_exp_and_mut_with_gene_alias.RData")
# 
# #Gene list
# inputData <-  read.delim("Datasets/example_data_for_test_function/gene_mut/data_mutations_mskcc_maf.txt", sep = "\t", stringsAsFactors = F)
# inputData <- unique(inputData$Hugo_Symbol)
# #GeneList <- unlist(strsplit(input$genelist, split = "\n"))
# 
# a <- ShinyDeepDR_input_mut(inputData = inputData,data.format = "LIST",alias.table = "ccle_exp_and_mut_with_gene_alias.RData")

ShinyDeepDR_input_mut <- function(inputData,data.format = c("MATRIX", "MAF","LIST"),sample_names = "Sample1",
                                  selectID=NA,
                                  alias.table = "ccle_exp_and_mut_with_gene_alias.RData"){
  
  source("function/CheckGeneSymbol.R")
  load("data/ccle_exp_and_mut_index_missingVal.RData")
  
  library(methods)
  
  #created S4 object
  setClass(Class = "input_mut",
           representation = representation(inputData = "data.frame",
                                           inputData.type = "character",
                                           sample_ID = "character",
                                           genelist = "character",
                                           genelist.update = "data.frame",
                                           mutation_table = "data.frame"))
  #txt file
  if (toupper(data.format) == "MATRIX") {
    
    colnames(inputData)[1] <- "Gene"
    if (!is.na(selectID)) {
      inputData <- inputData[,c("Gene", selectID)]
    }
    
    if(class(inputData$Gene) == "character"){
      
      genelist.update = ShinyDeepDR_CheckGeneSymbol(inputData = inputData,data.type = "mut",alias.table = alias.table)
      mutation_table.out <- merge(genelist.update, inputData, by.x = "gene.symbol", by.y = "Gene", all.x = T, sort = F)
      mutation_table.out <- mutation_table.out[!is.na(mutation_table.out$gene.symbol.update), -1]
      mutation_table.out <- merge(ccle.mut.index.missingVal, mutation_table.out, by.x = "Gene", by.y = "gene.symbol.update", all.x = T, sort = F)
      mutation_table.out <- mutation_table.out[match(ccle.mut.index.missingVal$Gene,mutation_table.out$Gene ),]
      
      identical(ccle.mut.index.missingVal$Gene, mutation_table.out$Gene)
      
      idx <- which(is.na(rowSums(mutation_table.out[,-1])))
      for (i in idx) {
        
        id.sample <- which(is.na(mutation_table.out[i,]))
        mutation_table.out[i, id.sample] <- mutation_table.out$missing_val[i]
        
      }
      
      mut.out <- new("input_mut",
                     inputData = inputData,
                     inputData.type = toupper(data.format),
                     sample_ID = colnames(inputData)[-1],
                     genelist.update = genelist.update,
                     genelist = inputData$Gene,
                     mutation_table = mutation_table.out[,-2]
      )
      
      
    }else{
      
      cat("Data format error! First column is gene symbols. \n")
    }
    
  }
  
  #MAF file
  if (toupper(data.format) == "MAF" ) {
    if (!is.na(selectID)) {
      inputData <- inputData[which(inputData$Tumor_Sample_Barcode %in% selectID),]
    }
    
    if( sum(toupper(colnames(inputData)) %in% toupper(c("Hugo_Symbol" ,"Variant_Classification","Tumor_Sample_Barcode"  ))) == 3){
      tb <- inputData[, toupper(c("Hugo_Symbol" ,"Variant_Classification","Tumor_Sample_Barcode"  ))]
      variant_types <- toupper(c("Missense_Mutation", "Nonsense_Mutation",  "Frame_Shift_Ins", "Frame_Shift_Del"))
      tb <- tb[which(toupper( tb$VARIANT_CLASSIFICATION) %in% variant_types),]
      tb$HUGO_SYMBOL <- toupper(tb$HUGO_SYMBOL)
      
      #Generate mutation table
      mutation_table <- data.frame(matrix(data = 0, nrow = length(unique(tb$HUGO_SYMBOL)),
                                          ncol = length(unique(tb$TUMOR_SAMPLE_BARCODE))+1), 
                                   stringsAsFactors = F)
      
      
      colnames(mutation_table) <- c("Gene", unique(tb$TUMOR_SAMPLE_BARCODE))
      rownames(mutation_table)  <- sort(unique(tb$HUGO_SYMBOL), decreasing = F)
      mutation_table$Gene       <- rownames(mutation_table)
      
      for (i in 1:nrow(mutation_table)) {
        tb1 <- tb[which(tb$HUGO_SYMBOL == mutation_table$Gene[i] ), ]
        mutation_table[i,which(colnames(mutation_table) %in% unique(tb1$TUMOR_SAMPLE_BARCODE))] <- 1
        
        #Part of gene have deplicate variant types in one sample
        #  if (nrow(tb1) != sum(colnames(mutation_table) %in% unique(tb1$Tumor_Sample_Barcode))) {
        #    cat(mutation_table$Gene[i], "idx:", i ,"\ Samples: ", sum(colnames(mutation_table) %in% unique(tb1$Tumor_Sample_Barcode)), "\n")
        #  }
        # }
      }
      
      genelist.update = ShinyDeepDR_CheckGeneSymbol(inputData = mutation_table,data.type = "mut",alias.table = alias.table)
      
      mutation_table.out <- merge(genelist.update, mutation_table, by.x = "gene.symbol", by.y="Gene",all.y = T ,sort = F)
      mutation_table.out <-  mutation_table.out[!is.na(mutation_table.out$gene.symbol.update),-1]
      colnames(mutation_table.out)[1] <- "Gene"
      
      mutation_table.out <- merge(ccle.mut.index.missingVal, mutation_table.out, by="Gene", all.x = T, sort = F)
      idx <- rowSums(mutation_table.out[,-1])
      
      for (i in which(is.na(idx))) {
        mutation_table.out[i, -c(1,2)] <- mutation_table.out$missing_val[i]
        
      }
      
      mutation_table.out <- mutation_table.out[match(ccle.mut.index.missingVal$Gene,mutation_table.out$Gene ),]
      identical(ccle.mut.index.missingVal$Gene, mutation_table.out$Gene)
      
      mut.out <- new("input_mut",
                     inputData = inputData,
                     inputData.type = toupper(data.format),
                     sample_ID = unique(tb$TUMOR_SAMPLE_BARCODE),
                     genelist.update = genelist.update,
                     genelist = unique(tb$HUGO_SYMBOL),
                     mutation_table = mutation_table.out[,-2]
      )
      
      
    }else{
      
      idx <- which( toupper(c("Hugo_Symbol" ,"Variant_Classification","Tumor_Sample_Barcode"  ) )%in% toupper(colnames(inputData) ))
      
      cat( paste(toupper(c("Hugo_Symbol" ,"Variant_Classification","Tumor_Sample_Barcode"))[-idx], collapse  = " & "),
           "column missing!  \n")
    }
  }
  
  
  #Gene list
  if (toupper(data.format) == "LIST" ) {
    if(class(inputData) == "character"){
      mutation_table  <- data.frame(Gene = inputData, stringsAsFactors = F)
      genelist.update = ShinyDeepDR_CheckGeneSymbol(inputData =mutation_table ,data.type = "mut",alias.table = alias.table)
      
      mutation_table <- merge(mutation_table, genelist.update, by.x= "Gene", by.y = "gene.symbol",all.y = T, sort = F)
      mutation_table <- mutation_table[!is.na(mutation_table$gene.symbol.update),]
      mutation_table.out <- ccle.mut.index.missingVal
      mutation_table.out$Sample1 <- 0
      mutation_table.out$Sample1[which(mutation_table.out$Gene %in% mutation_table$gene.symbol.update)] <- 1
      mutation_table.out <- mutation_table.out[match(ccle.mut.index.missingVal$Gene,mutation_table.out$Gene ),]
      identical(ccle.mut.index.missingVal$Gene, mutation_table.out$Gene)
      colnames(mutation_table.out)[which(colnames(mutation_table.out) == "Sample1")] <- sample_names
      
      mut.out <- new("input_mut",
                     inputData.type = toupper(data.format),
                     sample_ID = sample_names,
                     genelist.update = genelist.update,
                     genelist = unique(inputData),
                     mutation_table = mutation_table.out[,-2]
      )
      
    }else{
      
      cat("Data format error!  \n")
    }
    
    
  }
  
  
  return(mut.out)
  
}




