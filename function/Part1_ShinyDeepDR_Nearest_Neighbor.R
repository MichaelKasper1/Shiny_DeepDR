#Part1: Nearest neighbors functions-expression
#Date: 2022.11.29
# setwd("~/Desktop/shinyDeepDR_work/")
# 
# #CCLE
# load("Datasets/CCLE_22Q2/ccle_22Q2_expression_cells_1406.RData")
# #EXP
# input.data <- ccle.exp.22Q2[,c(1,10:13)]
# data.select ="CCLE"
# select.sample = colnames(input.data)[3]
# data.source = "data/Part1_Nearest_Neighbors/gdsc_and_ccle_704_overlapped_cells_exp_tpm.RDS"
# drug.source.exp = "data/GDSC_paper/gdsc_and_ccle_704_overlapped_cells_ic50.RDS"
# drug.source.mut = "data/GDSC_paper/gdsc_265drug_info.RDS"
# drug.info = "data/GDSC_paper/gdsc_265drug_info.RDS"
# cell.info = "data/GDSC_paper/gdsc_and_ccle_704_overlapped_cells_info.RDS"
# 
# #TCGA
# load("Datasets/TCGA/tcga_expression_samples_9059.RData")
# #EXP
# input.data <- tcga.exp[,c(1,10:13)]
# data.select ="TCGA"
# select.sample = colnames(input.data)[4]
# data.source = "Datasets/Part1_Nearest_Neighbors/tcga_expression_samples_9059.RData"
# #drug.source = "Datasets/GDSC_paper/gdsc_and_ccle_704_overlapped_cells_ic50.RDS"
# #drug.info = "Datasets/GDSC_paper/gdsc_265drug_info.RDS"
# cell.info = "Datasets/GDSC_paper/gdsc_990_cell_info.RDS"

  
ShinyDeepDR_NearestNeighbor_cor <- function(exp.data = NULL,
                                            mut.data = NULL,
                                            select.sample = NA,
                                            data.select = c("TCGA", "CCLE","PDXE"),
                                            cell.info,
                                            mutation_similarity_method= c("Fisher","Jaccard"),
                                            data.source, #CCLE: tpm (exp) and mut (exp) of 704 overlapped cells./ TCGA: exp and mut of 9059
                                            drug.ic50, #RDS format (GDSC)
                                            drug.info #RDS format (GDSC)
                                            ){
  
  library('dplyr')
  library('tibble')
  library('methods')
  library('WGCNA')
  source("function/Part1_shinyDeepDR_Mut_Similarity.R")
  
  #created S4 object
  setClass(Class = "nearest_neighbors",
           representation = representation(inputData.type = "character",
                                           data.select = "character",
                                           cor_table ="data.frame",
                                           top3.cor.drug = "data.frame",
                                           mut_similarity = "data.frame",
                                           top3.fisher.drug = "data.frame",
                                           top3.jaccard.drug = "data.frame"))
  
 
  
  if(is.na(select.sample)){
    if (!is.null(exp.data)) {
      exp.data = exp.data[,c(1,2)]  
      rownames(exp.data) <- toupper(exp.data$Gene)
      select.sample <- colnames(exp.data)[2]
    }
    
    if(!is.null(mut.data)){
      mut.data = mut.data[,c(1,2)]  
      rownames(mut.data) <- toupper(mut.data$Gene)
      select.sample <- colnames(mut.data)[2]
     
    }
    
    # if( colnames(exp.data)[2] == colnames(mut.data)[2]){
    #   select.sample = colnames(exp.data)[2]
    # }else{
    #   warning("Sample names do not match!!")
    # }
    # 
  
  }else{
    if (!is.null(exp.data)) {
      exp.data = exp.data[,c("Gene", select.sample)]  
      rownames(exp.data) <- toupper(exp.data$Gene)
    }
    
    if(!is.null(mut.data)){
      mut.data = mut.data[,c("Gene", select.sample)]  
      rownames(mut.data) <- toupper(mut.data$Gene)
    }
    
  }
  
  
  #Cell lines
  if (toupper(data.select) == "CCLE") {
    #Cell lines information
    cell.info <- cell.info
    #Drug information
    drug.info <- drug.info
    
    #GDSC IC50 table with 704 overlapped cell lines (overlap with CCLE exp and mut)
    drug.ic50 <- drug.ic50
    rownames(drug.ic50) <- drug.ic50$DepMapID
    
    if (!is.null(exp.data)) {
      
      #CCLE ---> overlapped cell lines : 704
      compaired.table = data.source
      rownames(compaired.table) <- toupper(compaired.table$Gene)
      compaired.table <- compaired.table[,-1]
      
      if(sum(rownames(exp.data) == rownames(compaired.table)) != 15363){
        exp.data <- exp.data[match(rownames(compaired.table), rownames(exp.data)),]
      }
      
      #calculate correlation
      # Pearson Correlation
      cor_table <- as.data.frame(round(t(WGCNA::cor(x = exp.data[,2],y = compaired.table, use = "all.obs",
                                                    nThreads = 16, method = "pearson")),digits = 3)) #WGCNA::cor: calculate faster
      cor_table <- cor_table %>% add_column(rownames(cor_table),.before = "V1") 
      colnames(cor_table) <-  c("DepMapID","Corr")
      cor_table <- cor_table[order(cor_table$Corr, decreasing = T), ]
      
      #add cell info
      cor_table_cell <- merge(x = cor_table,
                              y = cell.info,
                              by.x = "DepMapID",
                              by.y="DepMapID",
                              all.x = T,
                              sort = F)
      
      #add drug info
      top3cell.drug <- drug.ic50[which(drug.ic50$DepMapID %in% cor_table$DepMapID[1:3]),]
      rownames(top3cell.drug) <- top3cell.drug$DepMapID
      top3cell.drug.t <- as.data.frame(t(top3cell.drug[,-c(1,2)]))
      top3cell.drug.t <- top3cell.drug.t %>% add_column(DrugNames = rownames(top3cell.drug.t), .before = colnames(top3cell.drug.t)[1])
      top3cell.drug.t <- top3cell.drug.t[,c("DrugNames", cor_table_cell$DepMapID[1:3])]
      colnames(top3cell.drug.t)[2:4] <- paste0( cor_table_cell$cell_line_name[1:3], " (", cor_table_cell$DepMapID[1:3],")")
      top3cell.drug.out <- merge(x = top3cell.drug.t,
                                 y = drug.info,
                                 by.x = "DrugNames",
                                 by.y="colname",
                                 all.x = T,
                                 sort = F)
      
      final <- new("nearest_neighbors",
                   inputData.type = "exp",
                   data.select = data.select,
                   cor_table = cor_table_cell,
                   top3.cor.drug =  top3cell.drug.out)
      
      
    }
      
    
    if (!is.null(mut.data)) {
      
      #CCLE ---> overlapped cell lines : 704
      compaired.table = data.source
      compaired.table$Gene <- toupper(compaired.table$Gene)
      rownames(compaired.table) <- compaired.table$Gene
      
      if(sum(rownames(mut.data) == rownames(compaired.table)) != 18281){
        mut.data <- mut.data[match(toupper(rownames(compaired.table)), toupper(rownames(mut.data))),]
      }
      
      #calculate similarity
      if (tolower(mutation_similarity_method) == "fisher") {
        
      #Fisher exact test
      fisher_table <- ShinyDeepDR_Mut_Similarity(input.data = mut.data,
                                                  select.sample = select.sample,method = "fisher",
                                                  compaired.table = compaired.table)
      fisher_table <- fisher_table[order(fisher_table$value, decreasing = F),]
      
      #add cell info
      fisher_table_cell <- merge(x = fisher_table,
                                 y = cell.info,
                                 by.x = "DepMapID",
                                 by.y="DepMapID",
                                 all.x = T,
                                 sort = F)
      #add drug info
      top3cell.drug.fisher <- drug.ic50[which(drug.ic50$DepMapID %in% fisher_table_cell$DepMapID[1:3]),]
      rownames(top3cell.drug.fisher ) <- top3cell.drug.fisher$DepMapID
      top3cell.drug.fisher.t <- as.data.frame(t(top3cell.drug.fisher[,-c(1,2)]))
      top3cell.drug.fisher.t <- top3cell.drug.fisher.t %>% add_column(DrugNames = rownames(top3cell.drug.fisher.t), .before = colnames(top3cell.drug.fisher.t)[1])
      top3cell.drug.fisher.t <- top3cell.drug.fisher.t[,c("DrugNames", fisher_table_cell$DepMapID[1:3])]
      colnames(top3cell.drug.fisher.t)[2:4] <- paste0(fisher_table_cell$cell_line_name[1:3], " (", fisher_table_cell$DepMapID[1:3],")")
      top3cell.drug.fisher.out <- merge(x =top3cell.drug.fisher.t,
                                        y = drug.info,
                                        by.x = "DrugNames",
                                        by.y="colname",
                                        all.x = T,
                                        sort = F)
      
      
      final <- new("nearest_neighbors",
                   inputData.type = "mut",
                   data.select = data.select,
                   mut_similarity = fisher_table_cell,
                   top3.fisher.drug = top3cell.drug.fisher.out
      )
      
      }
      
      if (tolower(mutation_similarity_method) == "jaccard") {
        #Jaccard similarity
        Jaccard_table <- ShinyDeepDR_Mut_Similarity(input.data = mut.data,
                                                   select.sample = select.sample,
                                                   method = "jaccard",
                                                   compaired.table = compaired.table)
          
        #add cell info
        jaccard_table_cell <- merge(x = Jaccard_table,
                                   y = cell.info,
                                   by.x = "DepMapID",
                                   by.y="DepMapID",
                                   all.x = T,
                                   sort = F)
        jaccard_table_cell <- jaccard_table_cell[order(jaccard_table_cell$value, decreasing = T),]
        jaccard_table_cell$value <- round(jaccard_table_cell$value, digits = 3)
        
        
        #add drug info
        top3cell.drug.jaccard <- drug.ic50[which(drug.ic50$DepMapID %in% jaccard_table_cell$DepMapID[1:3]),]
        rownames(top3cell.drug.jaccard ) <- top3cell.drug.jaccard$DepMapID
        top3cell.drug.jaccard.t <- as.data.frame(t(top3cell.drug.jaccard[,-c(1,2)]))
        top3cell.drug.jaccard.t <- top3cell.drug.jaccard.t %>% add_column(DrugNames = rownames(top3cell.drug.jaccard.t), .before = colnames(top3cell.drug.jaccard.t)[1])
        top3cell.drug.jaccard.t <- top3cell.drug.jaccard.t[,c("DrugNames", jaccard_table_cell$DepMapID[1:3])]
        colnames(top3cell.drug.jaccard.t)[2:4] <- paste0(jaccard_table_cell$cell_line_name[1:3], " (", jaccard_table_cell$DepMapID[1:3],")")
        top3cell.drug.jaccard.out <- merge(x =top3cell.drug.jaccard.t,
                                          y = drug.info,
                                          by.x = "DrugNames",
                                          by.y="colname",
                                          all.x = T,
                                          sort = F)
        
        
        final <- new("nearest_neighbors",
                     inputData.type = "mut",
                     data.select = data.select, 
                     mut_similarity = jaccard_table_cell,
                     top3.jaccard.drug = top3cell.drug.jaccard.out
        )
        
      }
      
      
    }

  }
  
  ######################################
  #TCGA tumor
  if (toupper(data.select) == "TCGA") {
    #TCGA patients information
    cell.info <- cell.info
    
    if (!is.null(exp.data)) {

      #Input data type: Gene expression or both
      if(class(data.source) == "character"){
        load(data.source)
      compaired.table <- tcga.exp
      rownames(compaired.table) <- compaired.table$Gene

    
      }else{
        compaired.table = data.source
        rownames(compaired.table) <- compaired.table$Gene
      }
      
      if (colnames(compaired.table)[2] == "missing_val") {
        compaired.table <- compaired.table[,-2]
      }
      
      
      if(sum(rownames(exp.data) == rownames(compaired.table)) != 15363){
        exp.data <- exp.data[match(rownames(compaired.table), rownames(exp.data)),]
      }
      
    }
    
    if (!is.null(mut.data)) {
        
        #Input data type: Gene expression or both
        if(class(data.source) == "character"){
          load(data.source)
          compaired.table <- tcga.mut
          rownames(compaired.table) <- compaired.table$Gene
          
          
        }else{
          
          compaired.table = data.source
          rownames(compaired.table) <- compaired.table$Gene
        }
        
        if (colnames(compaired.table)[2] == "missing_val") {
          compaired.table <- compaired.table[,-2]
        }
      
      if(sum(rownames(mut.data) == rownames(compaired.table)) != 18281){
        mut.data <- mut.data[match(rownames(compaired.table), rownames(mut.data)),]
      }
      
    }
    
    
    #calculate correlation
    if (!is.null(exp.data)) {
      # Pearson Correlation
      cor_table <- as.data.frame(round(t(WGCNA::cor(x = exp.data[,2],y = compaired.table[,-1], use = "all.obs", nThreads = 16,method = "pearson")),digits = 3)) #WGCNA::cor: calculate faster
      cor_table <- cor_table %>% add_column(rownames(cor_table),.before = "V1") 
      colnames(cor_table) <-  c("PatientID","Corr")
      cor_table$PatientID <- substr(cor_table$PatientID,1,12)
      cor_table <- cor_table[order(cor_table$Corr, decreasing = T), ]
      
      #add patient clinical data
      cor_table_cell <- merge(x = cor_table,
                              y = cell.info,
                              by.x = "PatientID",
                              by.y="bcr_patient_barcode",
                              all.x = T,
                              sort = F)
      
      final <- new("nearest_neighbors",
                   inputData.type = "exp",
                   data.select = data.select,
                   cor_table = cor_table_cell)
 
    }
    
    if (!is.null(mut.data)) {
      
      #calculate similarity
      if (tolower(mutation_similarity_method) == "fisher") {
        
      #Fisher exact test
      fisher_table <- ShinyDeepDR_Mut_Similarity(input.data = mut.data,
                                                   select.sample = select.sample,method = "fisher",
                                                   compaired.table = compaired.table)
      fisher_table$DepMapID <- substr(fisher_table$DepMapID,1,12)
      colnames(fisher_table)[1] <-"PatientID"
      fisher_table <- fisher_table[order(fisher_table$value, decreasing = F),]
      
      #add patient clinical data
      fisher_table_cell <- merge(x = fisher_table,
                                 y = cell.info,
                                 by.x = "PatientID",
                                 by.y="bcr_patient_barcode",
                                 all.x = T,
                                 sort = F)
      
      
      final <- new("nearest_neighbors",
                   inputData.type = "mut",
                   data.select = data.select,
                   mut_similarity = fisher_table_cell)
    
      
      }
      
      if (tolower(mutation_similarity_method) == "jaccard") {
        #Jaccard similarity
        Jaccard_table <- ShinyDeepDR_Mut_Similarity(input.data = mut.data,
                                                    select.sample = select.sample,
                                                    method = "jaccard",
                                                    compaired.table = compaired.table)
        Jaccard_table$DepMapID <- substr(Jaccard_table$DepMapID,1,12)
        colnames(Jaccard_table)[1] <-"PatientID"
        
        #add patient clinical data
        jaccard_table_cell <- merge(x = Jaccard_table,
                                    y = cell.info,
                                    by.x = "PatientID",
                                    by.y="bcr_patient_barcode",
                                    all.x = T,
                                    sort = F)
        jaccard_table_cell <- jaccard_table_cell[order(jaccard_table_cell$value, decreasing = T),]
        jaccard_table_cell$value <- round(jaccard_table_cell$value, digits = 3)
        
        final <- new("nearest_neighbors",
                     inputData.type = "mut",
                     data.select = data.select,
                     mut_similarity = jaccard_table_cell)
        
        }
      
      
    
    }

    
    }
  
  
  
  return(final)
  
  
}
