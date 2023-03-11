#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library("rintrojs")
library("shiny")
library("shinyjs")
library("shinyBS")
library("shinydashboardPlus")
library("shinycssloaders")
library("shinydashboard")
library("shinyWidgets")
library("plotly")
library("DT")
library("VennDiagram")
library("rcdk")
library("tidyr")
library("stringr")
library("dplyr")
#library("sqldf")
library("igraph")
library("data.table")
library("visNetwork")
library("dplyr")
library("igraph")
library("tidyverse")
library("tensorflow")
library("keras")
library("WGCNA")
library('tibble')
library('methods')

#library(BiocManager)
#options(repos = BiocManager::repositories())

#library('reticulate')
# virtualenv_create("r-reticulate", python = "C:\\Users\\tapsya\\AppData\\Local\\Programs\\Python\\Python39")


options(shiny.maxRequestSize = 10*1024^2)

#setwd("~/OneDrive - University of Pittsburgh/shinyDeepDR/DeepDRShiny/")

source("function/ShinyDeepDR_input_exp.R")
source("function/ShinyDeepDR_input_mut.R")
source("function/ShinyDeepDR_Check_Sample_ID_exp_and_mut.R")
source("function/CheckGeneSymbol.R")
source("function/gene_exp_Convert2TPM.R")
source("function/Predictions_Functions.R")
source("function/Part1_shinyDeepDR_Mut_Similarity.R")
source("function/Part1_ShinyDeepDR_Nearest_Neighbor.R")
source("function/Part2_Prediction_output_function.R")
source("function/ShinyDeepDR_visNetwork_nn.R")

#load data
load("data/drug_names.RData")
load("data/ccle_exp_and_mut_index_missingVal.RData")
load("data/ccle_exp_and_mut_with_gene_alias.RData")

#GDSC cell info and drug info
cell.info = readRDS("data/GDSC_paper/gdsc_and_ccle_704_overlapped_cells_info.RDS")

data.source.exp = readRDS("data/GDSC_paper/gdsc_and_ccle_704_overlapped_cells_exp_tpm.RDS")
data.source.mut = readRDS("data/GDSC_paper/gdsc_and_ccle_704_overlapped_cells_mut.RDS")
drug.ic50 = readRDS("data/GDSC_paper/gdsc_and_ccle_704_overlapped_cells_ic50.RDS")
drug.info = readRDS("data/GDSC_paper/gdsc_265drug_info_update.RDS")


#TCGA data
load("data/TCGA/tcga_expression_samples_9059.RData") #TPM values form Chris
load("data/TCGA/tcga_mutation_samples_9059.RData")
tcga_ic50Pred = readRDS("data/TCGA/tcga_predict_ic50_9059samples.RDS")
tcga_clinicalData = readRDS('data/TCGA/TCGA_PanCan_9059_Clinical_data.RDS')

#load CCLE data
load("data/CCLE_22Q2/ccle_22Q2_expression_cells_1406.RData")
load("data/CCLE_22Q2/ccle_22Q2_mutation_cells_1771.RData")
ccle_ic50_Pred <-
  readRDS("data/CCLE_22Q2/CCLE_704_cells_exp_ic50_prediction.RDS")


#load models
model <- load_model_hdf5("Predictions/model_final.h5")
model2 <-
  load_model_hdf5("Predictions/model_final_exp.h5") #exp only model
model3 <-
  load_model_hdf5("Predictions/model_final_mut.h5") #mutation only model

#example data
ex.exp <-
  readRDS("data/example_data_for_test_function/gene_exp/gene_exp_tcga_tpm_1samples.RDS")
ex.mut <-
  readRDS("data/example_data_for_test_function/gene_mut/gene_mut_tcga_maf_1samples.RDS")

#Load help page and privacy
helpInfo <- read.delim("www/Figure_explain.csv", sep = ",")
#privacytag <- read.

#Load function
Pred_Gen <- function(mut_data, exp_data, model) {
  sample_final <- list(as.matrix(mut_data), as.matrix(exp_data))
  predictions <- predict(model, sample_final)
  predictions <- round(predictions, 4)
  return(data.frame(predictions))
}

Pred_Gen2 <- function(mut_data, model) {
  sample_final <- list(as.matrix(mut_data))
  predictions <- predict(model, sample_final)
  predictions <- round(predictions, 4)
  return(data.frame(predictions))
}

Pred_Gen3 <- function(exp_data, model) {
  sample_final <- list(as.matrix(exp_data))
  predictions <- predict(model, sample_final)
  predictions <- round(predictions, 4)
  return(data.frame(predictions))
}


# Define server logic
shinyServer(function(input, output, session) {
  ####Privacy policy:#####
  observeEvent(input$privacytag1, {
    shiny::showModal(shiny::modalDialog(
      size = "l",
      includeHTML("www/PrivacyPolicy.html"),
      easyClose = TRUE
    ))
  })
  
  
  #help tab content
  #help page for input data format
  observeEvent(input$help_input1
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab3.html"),
                   easyClose = TRUE
                 ))
               })
  
  observeEvent(input$help_input2
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab3.html"),
                   easyClose = TRUE
                 ))
               })
  
  observeEvent(input$helptab3p1
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab3.html"),
                   
                   easyClose = TRUE,
                 ))
               })
  
  observeEvent(input$helptab3p2
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab3.html"),
                   
                   easyClose = TRUE,
                 ))
               })
  
  observeEvent(input$helptab3p3
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab3.html"),
                   
                   easyClose = TRUE,
                 ))
               })
  
  observeEvent(input$helptab3p4
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab3.html"),
                   
                   easyClose = TRUE,
                 ))
               })
  
  observeEvent(input$helptab13p5
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab3.html"),
                   
                   easyClose = TRUE,
                 ))
               })
  
  #help page for Find Similar Samples (Find Sample)
  observeEvent(input$helptab1p1
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab1.html"),
                   
                   easyClose = TRUE,
                 ))
               })
  
  observeEvent(input$helptab1p2
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab1.html"),
                   
                   easyClose = TRUE,
                 ))
               })
  
  observeEvent(input$helptab1p3
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab1.html"),
                   
                   easyClose = TRUE,
                 ))
               })
  
  observeEvent(input$helptab1p4
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab1.html"),
                   
                   easyClose = TRUE,
                 ))
               })
  
  observeEvent(input$helptab1p5
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab1.html"),
                   
                   easyClose = TRUE,
                 ))
               })
  
  #help page for Predict Drug Response (Find Drug)
  observeEvent(input$helptab2p1
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab2.html"),
                   
                   easyClose = TRUE,
                 ))
               })
  
  observeEvent(input$helptab2p2
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab2.html"),
                   
                   easyClose = TRUE,
                 ))
               })
  
  observeEvent(input$helptab2p3
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab2.html"),
                   
                   easyClose = TRUE,
                 ))
               })
  
  observeEvent(input$helptab2p4
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab2.html"),
                   
                   easyClose = TRUE,
                 ))
               })
  
  observeEvent(input$helptab2p5
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab2.html"),
                   
                   easyClose = TRUE,
                 ))
               })
  
  
  
  #Question mark content
  #Figure2-1
  observeEvent(input$nnboth1, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[8]),
      
      easyClose = TRUE
    ))
  })
  
  #Figure2-2
  observeEvent(input$nnboth3, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[9]),
      
      easyClose = TRUE
    ))
  })
  
  #Figure2-1
  observeEvent(input$nnmut1, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[8]),
      
      easyClose = TRUE
    ))
  })
  
  #Figure2-2
  observeEvent(input$nnmut3, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[9]),
      
      easyClose = TRUE
    ))
  })
  
  #Figure2-1
  observeEvent(input$nnexp1, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[8]),
      
      easyClose = TRUE
    ))
  })
  
  #Figure2-2
  observeEvent(input$nnexp3, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[9]),
      
      easyClose = TRUE
    ))
  })
  
  
  #Predictions____________________
  
  observeEvent(input$pdboth1, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[1]),
      easyClose = TRUE
    ))
    
  })
  observeEvent(input$pdboth2, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[2]),
      
      easyClose = TRUE
    ))
  })
  observeEvent(input$pdboth3, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[3]),
      
      easyClose = TRUE
    ))
  })
  observeEvent(input$pdboth5, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[4]),
      
      easyClose = TRUE
    ))
  })
  observeEvent(input$pdboth6, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[5]),
      
      easyClose = TRUE
    ))
  })
  
  #NEW MODAL
  observeEvent(input$pdboth6sub, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[6]),
      
      easyClose = TRUE
    ))
  })
  observeEvent(input$pdboth7, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[7]),
      
      easyClose = TRUE
    ))
  })
  
  observeEvent(input$pdmut1, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[1]),
      
      easyClose = TRUE
    ))
  })
  observeEvent(input$pdmut2, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[2]),
      
      easyClose = TRUE
    ))
  })
  observeEvent(input$pdmut3, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[3]),
      
      easyClose = TRUE
    ))
  })
  observeEvent(input$pdmut5, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[4]),
      
      easyClose = TRUE
    ))
  })
  observeEvent(input$pdmut6, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[5]),
      
      easyClose = TRUE
    ))
  })
  
  #NEW MODAL
  observeEvent(input$pdmut6sub, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[6]),
      
      easyClose = TRUE
    ))
  })
  observeEvent(input$pdmut7, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[7]),
      
      easyClose = TRUE
    ))
  })
  
  observeEvent(input$pdexp1, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[1]),
      
      easyClose = TRUE
    ))
  })
  observeEvent(input$pdexp2, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[2]),
      
      easyClose = TRUE
    ))
  })
  observeEvent(input$pdexp3, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[3]),
      
      easyClose = TRUE
    ))
  })
  observeEvent(input$pdexp5, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[4]),
      
      easyClose = TRUE
    ))
  })
  observeEvent(input$pdexp6, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[5]),
      
      easyClose = TRUE
    ))
  })
  
  #NEW MODAL
  observeEvent(input$pdexp6sub, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[6]),
      
      easyClose = TRUE
    ))
  })
  observeEvent(input$pdexp7, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$html[7]),
      
      easyClose = TRUE
    ))
  })
  
  #####################################
  pred_data <- reactiveValues()
  pred_data$mut_data <- NA
  pred_data$exp_data <- NA
  
  pred_data$mut_data2 <- NA
  pred_data$exp_data2 <- NA
  
  pred_data$genelist <- NA
  pred_data$genelist2 <- NA
  #when data is fully uploaded, render the data table showing the predictions
  
  ################################################################################################
  #Module1:Nearest Neighbors:
  #####Upload both files#####
  # example or user data
  dat.source <- reactiveValues(source = NULL)
  
  #Example
  observeEvent(input$submitexample, {
    showModal(modalDialog(
      size = "m",
      HTML(
        "Analysis started. Please click “Dismiss” to view results. It make take a minute to run."
      ),
      easyClose = TRUE
    ))
    
    dat.source$source <- "example"
  })
  
  #Message for submit
  observeEvent(input$submitmodel, {
    showModal(modalDialog(
      size = "m",
      HTML(
        "Analysis started. Please click “Dismiss” to view results. It make take a minute to run."
      ),
      easyClose = TRUE
    ))
    dat.source$source <- "user"
  })
  
  #Read expression data
  observeEvent(input$exp_data, {
    exp <- input$exp_data
    exp_data <- read.delim(exp$datapath)
    
    if (class(exp_data[1, 1]) != "character") {
      # session$sendCustomMessage(type = "resetFileInputHandler", "exp_data")
      showModal(modalDialog(
        size = "m",
        HTML(
          "Please ensure that the expression data is the text file with tab-delimited foramt. <br>
        In addition, the format of file is
        the column being the samples and rows being the gene symbols. The first column should be the gene symbols."
        ),
        easyClose = TRUE
      ))
      
    } else if (length(colnames(exp_data)[-1]) == 1) {
      colnames(exp_data)[1] <- "Gene"
      exp_data$Gene <- toupper(exp_data$Gene)
      pred_data$exp_data = exp_data
      
      showModal(modalDialog(
        size = "m",
        HTML("The format of expression data is correct!"),
        easyClose = TRUE
      ))
      
    } else{
      showModal(modalDialog(
        size = "m",
        HTML(
          "The format of expression data is correct, but we only accept one sample!!!! "
        ),
        easyClose = TRUE
      ))
      
    }
    
    
  })
  
  #Mutation data by gene list
  observeEvent(input$gene_list, {
    genelist = toupper(unlist(strsplit(isolate(
      input$gene_list
    ), split = "\n")))
    genelist = genelist[which(genelist != "")]

    if (length(genelist) != 0) {
      mutation_table  <- data.frame(Gene = genelist, stringsAsFactors = F)
      genelist.update = ShinyDeepDR_CheckGeneSymbol(
        inputData = mutation_table ,
        data.type = "mut",
        alias.table = "data/ccle_exp_and_mut_with_gene_alias.RData"
      )
      
      if (sum(is.na(genelist.update$gene.symbol.update)) == 0) {
        pred_data$genelist <- genelist.update$gene.symbol.update
        
        showModal(modalDialog(
          size = "m",
          HTML("All gene symbols are valid."),
          easyClose = TRUE
        ))
        
      } else{
        pred_data$genelist <- genelist.update$gene.symbol.update[which(!is.na( genelist.update$gene.symbol.update))]
        showModal(modalDialog(size = "m",
                              HTML(
                                paste(
                                  "Invalid gene symbols : ",
                                  paste(c(genelist.update$gene.symbol[which(is.na(genelist.update$gene.symbol.update))]), collapse = "; "),
                                  "<br>",
                                  "These genes will be excluded from the analysis."
                                )
                              ),
                              easyClose = TRUE))
      }
      
    }
    
  })
  
  
  #Read mutation data (MAF files)
  observeEvent(input$mut_data, {
    mut <- input$mut_data
    mut_data <- read.delim(mut$datapath)
    colnames(mut_data) <- toupper(colnames(mut_data))
    
    if (sum(colnames(mut_data) %in% toupper(
      c(
        "Hugo_Symbol" ,
        "Variant_Classification",
        "Tumor_Sample_Barcode"
      )
    )) == 3) {
      mut_data = mut_data[, toupper(c(
        "Hugo_Symbol" ,
        "Variant_Classification",
        "Tumor_Sample_Barcode"
      ))]
      mut_data$HUGO_SYMBOL <- toupper(mut_data$HUGO_SYMBOL)
      pred_data$mut_data = mut_data
      session$sendCustomMessage(type = "resetFileInputHandler", "mut_data")
      
      if (length(unique(mut_data$TUMOR_SAMPLE_BARCODE)) != 1) {
        showModal(modalDialog(
          size = "m",
          HTML("Multiple samples!! <br> We only accept one sample!!"),
          easyClose = TRUE
          
        ))
        
      } else{
        showModal(modalDialog(
          size = "m",
          HTML(
            "If you already pasted a gene list, you do not need to submit the MAF file. Otherwise, the MAF file will take precedence the list."
          ),
          easyClose = TRUE
        ))
        
        
      }
      
      
    } else{
      session$sendCustomMessage(type = "resetFileInputHandler", "mut_data")
      showModal(modalDialog(
        size = "m",
        HTML(
          "Wrong format.<br>
        Please ensure that the mutation data is a MAF file or text file with Hugo_Symbol, Variant_Classification and
             Tumor_Sample_Barcode columns."
        ),
        easyClose = TRUE
      ))
      
      
    }
    
  })
  
  
  #table
  table.nn.both <-  reactive({
    if (is.null(dat.source$source))
      return()
    
    if (dat.source$source ==  "example") {
      exp_data =  ex.exp
      mut_data =  ex.mut
      genelist = NA
      
    }
    
    if (dat.source$source ==  "user") {
      exp_data =  pred_data$exp_data
      mut_data =  pred_data$mut_data
      
      if(!is.data.frame(mut_data)){
        genelist = pred_data$genelist
        
        print("Use genelist")
      }
      
      if(is.data.frame(mut_data)){
        genelist = NA
        print("Use mutation")
      }

    }
    
    data.values <- input$exp_type
    log_transform <- input$log_transform
    
    if (!is.data.frame(exp_data)) {
      showModal(modalDialog(
        size = "m",
        HTML("Please upload expression file"),
        easyClose = TRUE
      ))
      return()
      
    }
    
    if (!is.data.frame(mut_data) &&  is.na(genelist)) {
      showModal(modalDialog(
        size = "m",
        HTML("Please upload mutation file or paste gene list!"),
        easyClose = TRUE
      ))
      return()
    }
    
    
    if (is.data.frame(mut_data)) {
      #Check Sample names
      if (sum(colnames(exp_data)[2:ncol(exp_data)] != unique(mut_data$TUMOR_SAMPLE_BARCODE)) >
          0) {
        showModal(modalDialog(
          size = "m",
          HTML(
            "Sample name is not matched between expresstion and mutation data"
          ),
          easyClose = TRUE
          
        ))
      }
    }
    
    
    #Expression data
    exp_out <- ShinyDeepDR_input_exp(
      inputData = exp_data,
      data.values = data.values,
      selectID = NA,
      log2_transform = log_transform ,
      alias.table = "data/ccle_exp_and_mut_with_gene_alias.RData"
    )
    

    #Mutation data
    if (is.data.frame(mut_data) && is.na(genelist)) {
      mut_out <- ShinyDeepDR_input_mut(
        inputData = mut_data,
        data.format = "MAF",
        alias.table = "data/ccle_exp_and_mut_with_gene_alias.RData"
      )
    } else{
      mut_out <- ShinyDeepDR_input_mut(
        inputData = genelist,
        data.format = "LIST",
        sample_names = exp_out@sample_ID[1],
        alias.table = "data/ccle_exp_and_mut_with_gene_alias.RData"
      )
    }
    
    
    #Calculate Nearest neighbors
    cell_type = input$cell_type1
    
    if (cell_type == "CCLE") {
      ccle_nn_exp <-
        ShinyDeepDR_NearestNeighbor_cor(
          exp.data =  exp_out@exp_final,
          mut.data = NULL,
          data.select = "CCLE",
          cell.info = cell.info,
          data.source = data.source.exp,
          drug.ic50 = drug.ic50,
          drug.info = drug.info
        )
      ccle_exp <- ccle_nn_exp@cor_table
      ccle_exp <-
        ccle_exp %>% add_column(order_exp = seq(1, nrow(ccle_exp), 1), .after = "Corr")
      
      
      #mutation
      ccle_nn_mut <- ShinyDeepDR_NearestNeighbor_cor(
        exp.data = NULL,
        mut.data =  mut_out@mutation_table,
        data.select = "CCLE",
        mutation_similarity_method = "jaccard",
        cell.info = cell.info,
        data.source = data.source.mut,
        drug.ic50 = drug.ic50,
        drug.info = drug.info
      )
      ccle_mut <- ccle_nn_mut@mut_similarity
      ccle_mut <-
        ccle_mut %>% add_column(order_mut = seq(1, nrow(ccle_mut), 1), .after = "value")
      
      out_table =  merge(ccle_exp[, c(1:3)], ccle_mut[, c(1:3)], by = "DepMapID" , sort = F)
      out_table$order <- out_table$order_exp + out_table$order_mut
      out_table <- out_table[order(out_table$order, decreasing = F), ]
      out_table.final <- out_table[, c("DepMapID" , "Corr", "value")]
      out_table.final <-
        merge(out_table.final,
              cell.info,
              by = "DepMapID" ,
              sort = F)
      
    }
    
    if (cell_type == "TCGA") {
      #expression
      tcga_nn_exp <-
        ShinyDeepDR_NearestNeighbor_cor(
          exp.data =  exp_out@exp_final,
          mut.data = NULL,
          data.select = "TCGA",
          cell.info = tcga_clinicalData,
          data.source = tcga.exp
        )
      tcga_exp <- tcga_nn_exp@cor_table
      tcga_exp <-
        tcga_exp %>% add_column(order_exp = seq(1, nrow(tcga_exp), 1), .after = "Corr")
      
      
      #mutation
      tcga_nn_mut <- ShinyDeepDR_NearestNeighbor_cor(
        exp.data = NULL,
        mut.data = mut_out@mutation_table,
        data.select = "TCGA",
        mutation_similarity_method = "jaccard",
        cell.info = tcga_clinicalData,
        data.source = tcga.mut
      )
      tcga_mut <- tcga_nn_mut@mut_similarity
      tcga_mut <-
        tcga_mut %>% add_column(order_mut = seq(1, nrow(tcga_mut), 1), .after = "value")
      
      out_table =  merge(tcga_exp[, c(1:3)], tcga_mut[, c(1:3)], by = "PatientID" , sort = F)
      out_table$order <- out_table$order_exp + out_table$order_mut
      out_table <- out_table[order(out_table$order, decreasing = F), ]
      out_table.final <- out_table[, c("PatientID", "Corr", "value")]
      out_table.final <-
        merge(
          out_table.final,
          tcga_clinicalData,
          by.x = "PatientID",
          by.y = "bcr_patient_barcode" ,
          sort = F
        )
    }
    
    
    out_table.final
    
    #print(head(out_table.final))
    
  })
  
  #####Nearest_Neighbors_calculation#########
  
  # Table1-1: nn.table.both ---> nearest neighbors calculation and need to add sample information
  output$nn.table.both <- renderDataTable({
    cell_type = input$cell_type1
    table.out = table.nn.both()
    
    if (is.null(table.out))
      return()
    
    if (toupper(cell_type) == "CCLE") {
      #change column name
      colnames(table.out) <- c(
        "DepMapID",
        "Correlation Exp",
        "Jaccard index Mut",
        "Cell line name",
        "stripped cell line name",
        #hide
        "CCLE Name",
        "Alias",
        "COSMIC ID",
        "Sample collection site",
        "Primary or metastasis",
        "Primary disease" ,
        "Subtype",
        "Lineage",
        "Lineage subtype",
        "Lineage sub subtype",
        "Lineage molecular subtype",
        "Cellosaurus NCIt disease"
      )
      table.final <- table.out[, c(
        "Cell line name",
        "Jaccard index Mut",
        "Correlation Exp",
        "Primary disease",
        "Subtype",
        "Primary or metastasis",
        "Sample collection site",
        "Lineage",
        "Lineage subtype",
        "Lineage sub subtype",
        "Lineage molecular subtype",
        "Cellosaurus NCIt disease",
        "CCLE Name",
        "Alias",
        "DepMapID",
        "COSMIC ID"
      )]
      
    }
    
    if (toupper(cell_type) == "TCGA") {
      #change column name
      table.final <- table.out[, c(
        "PatientID",
        "value",
        "Corr",
        "type",
        "age_at_initial_pathologic_diagnosis",
        "gender",
        "race",
        "ajcc_pathologic_tumor_stage",
        "clinical_stage",
        "histological_type",
        "histological_grade"
      )]
      colnames(table.final) <- c(
        "Patient ID",
        "Jaccard index Mut",
        "Correlation Exp",
        "Cancer type",
        "Age at diagnosis",
        "Gender",
        "Race",
        "AJCC stage",
        "Clinical stage",
        "Histological type",
        "Histological grade"
      )
      
    }
    
    DT::datatable({
      table.final
    },
    extensions = c('Scroller', 'FixedColumns'), rownames = FALSE,
    options = list(
      pageLength = 10,
      lengthMenu = c(10, 15, 20),
      scrollX = TRUE
      # extensions = c('Buttons', 'FixedColumns'), rownames = FALSE,
      # rowCallback = JS(
      #   "function(row, data) {",
      #   "for (i = 2; i < 4; i++) {",
      #   "if (data[i]>1000 | data[i]<1){",
      #   "$('td:eq('+i+')', row).html(data[i].toExponential(3));",
      #   "}",
      #   "}",
      #   "}"
      #   )
      
      #     dom = 'Bfrtip',
      #     buttons = c('excel', 'csv'),
      #     scrollX = TRUE
    ),
    selection = list(mode = "single"))
    
  })
  
  # Figure1-2: nn.tsne.both
  
  # Figure1-3: nn.network.both
  output$nn.network.both <- renderVisNetwork({
    table.out = table.nn.both()
    cell_type = input$cell_type1
    
    if (is.null(table.out))
      return()
    
    ndx <- input$nn.table.both_rows_current
    if (is.null(ndx)) {
      ndx = 1:10
    }
    
    
    
    p <- ShinyDeepDR_visNetwork_nn(
      inputTable = table.out[ndx,],
      target = "Query_sample",
      cell_type = cell_type,
      gdsc_ic50 = drug.ic50 ,
      tcga.predict50 = tcga_ic50Pred,
      num_drug = 10
    )
    
    p
  })
  
  # Figure1-4: nn.ic50.both
  
  ##############################################################################################################
  #Upload mutation only
  # example or user data
  dat.source2 <- reactiveValues(source = NULL)
  #Example
  observeEvent(input$submitexample2, {
    showModal(modalDialog(
      size = "m",
      HTML(
        "Analysis started. Please click “Dismiss” to view results. It make take a minute to run."
      ),
      easyClose = TRUE
    ))
    
    dat.source2$source <- "example"
  })
  
  #Message for submit
  observeEvent(input$submitmodel2, {
    showModal(modalDialog(
      size = "m",
      HTML(
        "Analysis started. Please click “Dismiss” to view results. It make take a minute to run."
      ),
      easyClose = TRUE
    ))
    
    dat.source2$source <- "user"
  })
  
  #Mutation data by gene list
  observeEvent(input$gene_list2, {
    genelist = toupper(unlist(strsplit(input$gene_list2, split = "\n")))
    genelist = genelist[which(genelist != "")]
    
    if (length(genelist) != 0) {
      mutation_table  <- data.frame(Gene = genelist, stringsAsFactors = F)
      genelist.update = ShinyDeepDR_CheckGeneSymbol(
        inputData = mutation_table ,
        data.type = "mut",
        alias.table = "data/ccle_exp_and_mut_with_gene_alias.RData"
      )
      
      if (sum(is.na(genelist.update$gene.symbol.update)) == 0) {
        pred_data$genelist2 <- genelist.update$gene.symbol.update
        showModal(modalDialog(
          size = "m",
          HTML("All gene symbols are valid."),
          easyClose = TRUE
        ))
        
      } else{
        pred_data$genelist2 <- genelist.update$gene.symbol.update[which(!is.na( genelist.update$gene.symbol.update))]
        showModal(modalDialog(size = "m",
                              HTML(
                                paste(
                                  "Invalid gene symbols : ",
                                  paste(c(genelist.update$gene.symbol[which(is.na(genelist.update$gene.symbol.update))]), collapse = "; "),
                                  "<br>",
                                  "These genes will be excluded from the analysis."
                                )
                              ),
                              easyClose = TRUE))
      }
      
      
    }
    
    
  })
  
  
  #Read mutation data (MAF files)
  observeEvent(input$mut_data2, {
    mut <- input$mut_data2
    mut_data <- read.delim(mut$datapath)
    colnames(mut_data) <- toupper(colnames(mut_data))
    
    if (sum(toupper(colnames(mut_data)) %in% toupper(
      c(
        "Hugo_Symbol" ,
        "Variant_Classification",
        "Tumor_Sample_Barcode"
      )
    )) == 3) {
      mut_data = mut_data[, toupper(c(
        "Hugo_Symbol" ,
        "Variant_Classification",
        "Tumor_Sample_Barcode"
      ))]
      mut_data$HUGO_SYMBOL <- toupper(mut_data$HUGO_SYMBOL)
      pred_data$mut_data2 = mut_data
      
      if (length(unique(mut_data$TUMOR_SAMPLE_BARCODE)) != 1) {
        showModal(modalDialog(
          size = "m",
          HTML("Multiple samples!! <br> We only accept one sample!!"),
          easyClose = TRUE
          
        ))
        
      } else{
        showModal(modalDialog(
          size = "m",
          HTML(
            "If you already pasted a gene list, you do not need to submit the MAF file. Otherwise, the MAF file will take precedence the list."
          ),
          easyClose = TRUE
        ))
        
        
      }
      
      
    } else{
      session$sendCustomMessage(type = "resetFileInputHandler", "mut_data")
      showModal(modalDialog(
        size = "m",
        HTML(
          "Wrong format.<br>
        Please ensure that the mutation data is a MAF file or text file with Hugo_Symbol, Variant_Classification and
             Tumor_Sample_Barcode columns."
        ),
        easyClose = TRUE
      ))
      
      
    }
    
  })
  
  
  #####Nearest_Neighbors_calculation#########
  
  table.nn.mut <- reactive({
    if (is.null(dat.source2$source))
      return()
    
    if (dat.source2$source ==  "example") {
      mut_data =  ex.mut
      genelist = NA
    }
    
    
    if (dat.source2$source ==  "user") {
      mut_data =  pred_data$mut_data2
      
      if(!is.data.frame(mut_data)){
        genelist = pred_data$genelist2
        
        print("Use genelist")
      }
      
      if(is.data.frame(mut_data)){
        genelist = NA
        print("Use mutation")
      }
      
    }
    
    #Mutation data
    if (is.data.frame(mut_data) && is.na(genelist)) {
      mut_out <- ShinyDeepDR_input_mut(
        inputData = mut_data,
        data.format = "MAF",
        alias.table = "data/ccle_exp_and_mut_with_gene_alias.RData"
      )
    } else{
      mut_out <- ShinyDeepDR_input_mut(
        inputData = genelist,
        data.format = "LIST",
        sample_names = "Sample1",
        alias.table = "data/ccle_exp_and_mut_with_gene_alias.RData"
      )
    }
    
    
    
    #Calculate Nearest neighbors
    cell_type2 = input$cell_type2
    
    # #Calculate Nearest neighbors
    if (cell_type2 == "CCLE") {
      #mutation
      ccle_nn_mut <-
        ShinyDeepDR_NearestNeighbor_cor(
          exp.data = NULL,
          mut.data =  mut_out@mutation_table,
          data.select = "CCLE",
          mutation_similarity_method = "jaccard",
          cell.info = cell.info,
          data.source = data.source.mut,
          drug.ic50 = drug.ic50,
          drug.info = drug.info
        )
      out_table <- ccle_nn_mut@mut_similarity
      
    }
    
    if (cell_type2 == "TCGA") {
      #mutation
      tcga_nn_mut <-
        ShinyDeepDR_NearestNeighbor_cor(
          exp.data = NULL,
          mut.data = mut_out@mutation_table,
          data.select = "TCGA",
          mutation_similarity_method = "jaccard",
          cell.info = tcga_clinicalData,
          data.source = tcga.mut
        )
      out_table <- tcga_nn_mut@mut_similarity
      
    }
    
    out_table
    
  })
  
  # Table1-1: nn.table.mut  ---> nearest neighbors calculation and need to add sample information
  output$nn.table.mut <- renderDataTable({
    table.out = table.nn.mut()
    cell_type = input$cell_type2
    
    if (is.null(table.out))
      return()
    
    if (toupper(cell_type) == "CCLE") {
      #change column name
      colnames(table.out) <- c(
        "DepMapID",
        "Jaccard index Mut",
        "Cell line name",
        "stripped cell line name",
        #hide
        "CCLE Name",
        "Alias",
        "COSMIC ID",
        "Sample collection site",
        "Primary or metastasis",
        "Primary disease" ,
        "Subtype",
        "Lineage",
        "Lineage subtype",
        "Lineage sub subtype",
        "Lineage molecular subtype",
        "Cellosaurus NCIt disease"
      )
      table.final <- table.out[, c(
        "Cell line name",
        "Jaccard index Mut",
        "Primary disease",
        "Subtype",
        "Primary or metastasis",
        "Sample collection site",
        "Lineage",
        "Lineage subtype",
        "Lineage sub subtype",
        "Lineage molecular subtype",
        "Cellosaurus NCIt disease",
        "CCLE Name",
        "Alias",
        "DepMapID",
        "COSMIC ID"
      )]
      
    }
    
    if (toupper(cell_type) == "TCGA") {
      #change column name
      table.final <- table.out[, c(
        "PatientID",
        "value",
        "type",
        "age_at_initial_pathologic_diagnosis",
        "gender",
        "race",
        "ajcc_pathologic_tumor_stage",
        "clinical_stage",
        "histological_type",
        "histological_grade"
      )]
      colnames(table.final) <- c(
        "Patient ID",
        "Jaccard index Mut",
        "Cancer type",
        "Age at diagnosis",
        "Gender",
        "Race",
        "AJCC stage",
        "Clinical stage",
        "Histological type",
        "Histological grade"
      )
      
    }
    
    DT::datatable({
      table.final
    },
    extensions = c('Buttons', 'FixedColumns'), rownames = FALSE, escape = FALSE,
    options = list(
      dom = 'Bfrtip',
      buttons = c('excel', 'csv'),
      scrollX = TRUE
      # rowCallback = JS(
      #   "function(row, data) {",
      #   "for (i = 2; i < 3; i++) {",
      #   #"if (data[i]>1000 | data[i]<1){",
      #   "$('td:eq('+i+')', row).html(data[i].toExponential(3));",
      #   #"}",
      #   "}",
      #   "}"
      # )
    ),
    selection = list(mode = "single"))
    
  })
  
  
  # Figure1-2: nn.tsne.mut
  # Figure1-3: nn.network.mut
  output$nn.network.mut <- renderVisNetwork({
    table.out = table.nn.mut()
    cell_type = input$cell_type2
    
    if (is.null(table.out))
      return()
    
    ndx <- input$nn.table.mut_rows_current
    if (is.null(ndx)) {
      ndx = 1:10
    }
    
    
    p <- ShinyDeepDR_visNetwork_nn(
      inputTable = table.out[ndx,],
      target = "Query_sample",
      cell_type = cell_type,
      gdsc_ic50 = drug.ic50 ,
      tcga.predict50 = tcga_ic50Pred,
      num_drug = 10
    )
    
    p
  })
  # Figure1-4: nn.ic50.mut
  
  
  ##############################################################################################################
  #Upload exp data only
  
  # example or user data
  dat.source3 <- reactiveValues(source = NULL)
  
  #Read expression data
  observeEvent(input$exp_data2, {
    exp <- input$exp_data2
    exp_data <- read.delim(exp$datapath)
    
    if (class(exp_data[1, 1]) != "character") {
      session$sendCustomMessage(type = "resetFileInputHandler", "exp_data2")
      showModal(modalDialog(
        size = "m",
        HTML(
          "Please ensure that the expression data is the text file with tab-delimited foramt. <br>
        In addition, the format of file is
        the column being the samples and rows being the gene symbols. The first column should be the gene symbols."
        ),
        easyClose = TRUE
      ))
    } else{
      colnames(exp_data)[1] <- "Gene"
      exp_data$Gene <- toupper(exp_data$Gene)
      pred_data$exp_data2 = exp_data
      showModal(modalDialog(
        size = "m",
        HTML("The format of expression data is correct!"),
        easyClose = TRUE
      ))
    }
    
    
  })
  
  #Message for submit
  observeEvent(input$submitmodel3, {
    showModal(modalDialog(
      size = "m",
      HTML(
        "Analysis started. Please click “Dismiss” to view results. It make take a minute to run."
      ),
      easyClose = TRUE
    ))
    dat.source3$source <- "user"
  })
  #Example
  observeEvent(input$submitexample3, {
    showModal(modalDialog(
      size = "m",
      HTML(
        "Analysis started. Please click “Dismiss” to view results. It make take a minute to run."
      ),
      easyClose = TRUE
    ))
    
    dat.source3$source <- "example"
  })
  
  #table
  table.nn.exp <-  reactive({
    if (is.null(dat.source3$source))
      return()
    
    if (dat.source3$source ==  "example") {
      exp_data =  ex.exp
      
    }
    
    if (dat.source3$source ==  "user") {
      exp_data =  pred_data$exp_data2
      
      if (!is.data.frame(exp_data)) {
        showModal(modalDialog(
          size = "m",
          HTML("Please upload expression file"),
          easyClose = TRUE
        ))
        return()
        
      }
    }
    data.values <- input$exp_type2
    log_transform <- input$log_transform2
    #Expression data
    exp_out <- ShinyDeepDR_input_exp(
      inputData = exp_data,
      data.values = data.values,
      selectID = NA,
      log2_transform = log_transform ,
      alias.table = "data/ccle_exp_and_mut_with_gene_alias.RData"
    )
    
    
    
    # #Calculate Nearest neighbors
    cell_type = input$cell_type3
    
    if (cell_type == "CCLE") {
      #expression
      ccle_nn_exp <-
        ShinyDeepDR_NearestNeighbor_cor(
          exp.data =  exp_out@exp_final,
          mut.data = NULL,
          data.select = "CCLE",
          cell.info = cell.info,
          data.source = data.source.exp,
          drug.ic50 = drug.ic50,
          drug.info = drug.info
        )
      out_table <- ccle_nn_exp@cor_table
      
    }
    
    if (cell_type == "TCGA") {
      #expression
      tcga_nn_exp <-
        ShinyDeepDR_NearestNeighbor_cor(
          exp.data =  exp_out@exp_final,
          mut.data = NULL,
          data.select = "TCGA",
          cell.info = tcga_clinicalData,
          data.source = tcga.exp
        )
      out_table <- tcga_nn_exp@cor_table
      
    }
    
    
    
    out_table
    
    
  })
  
  #####Nearest_Neighbors_calculation#########
  
  # Table1-1: nn.table.exp ---> nearest neighbors calculation and need to add sample information
  output$nn.table.exp <- renderDataTable({
    cell_type = input$cell_type3
    table.out = table.nn.exp()
    
    if (is.null(table.out))
      return()
    
    if (toupper(cell_type) == "CCLE") {
      #change column name
      colnames(table.out) <- c(
        "DepMapID",
        "Correlation Exp",
        "Cell line name",
        "stripped cell line name",
        #hide
        "CCLE Name",
        "Alias",
        "COSMIC ID",
        "Sample collection site",
        "Primary or metastasis",
        "Primary disease" ,
        "Subtype",
        "Lineage",
        "Lineage subtype",
        "Lineage sub subtype",
        "Lineage molecular subtype",
        "Cellosaurus NCIt disease"
      )
      
      table.final <- table.out[, c(
        "Cell line name",
        "Correlation Exp",
        "Primary disease",
        "Subtype",
        "Primary or metastasis",
        "Sample collection site",
        "Lineage",
        "Lineage subtype",
        "Lineage sub subtype",
        "Lineage molecular subtype",
        "Cellosaurus NCIt disease",
        "CCLE Name",
        "Alias",
        "DepMapID",
        "COSMIC ID"
      )]
      
    }
    
    if (toupper(cell_type) == "TCGA") {
      #change column name
      table.final <- table.out[, c(
        "PatientID",
        "Corr",
        "type",
        "age_at_initial_pathologic_diagnosis",
        "gender",
        "race",
        "ajcc_pathologic_tumor_stage",
        "clinical_stage",
        "histological_type",
        "histological_grade"
      )]
      colnames(table.final) <- c(
        "Patient ID",
        "Correlation Exp",
        "Cancer type",
        "Age at diagnosis",
        "Gender",
        "Race",
        "AJCC stage",
        "Clinical stage",
        "Histological type",
        "Histological grade"
      )
      
    }
    
    DT::datatable({
      table.final
      
    },
    extensions = c('Scroller', 'FixedColumns'), rownames = FALSE,
    options = list(
      pageLength = 10,
      lengthMenu = c(10, 15, 20),
      scrollX = TRUE
      # extensions = c('Buttons', 'FixedColumns'), rownames = FALSE,
      # options = list(
      #     dom = 'Bfrtip',
      #     buttons = c('excel', 'csv'),
      #     scrollX = TRUE
    ),
    selection = list(mode = "single"))
    
  })
  
  # Figure1-2: nn.tsne.exp
  
  # Figure1-3: nn.network.exp
  output$nn.network.exp <- renderVisNetwork({
    table.out = table.nn.exp()
    cell_type = input$cell_type3
    
    if (is.null(table.out))
      return()
    
    ndx <- input$nn.table.exp_rows_current
    if (is.null(ndx)) {
      ndx = 1:10
    }
    
    
    p <- ShinyDeepDR_visNetwork_nn(
      inputTable = table.out[ndx,],
      target = "Query_sample",
      cell_type = cell_type,
      gdsc_ic50 = drug.ic50 ,
      tcga.predict50 = tcga_ic50Pred,
      num_drug = 10
    )
    
    p
  })
  
  # Figure1-4: nn.ic50.exp
  
  #Show description for each plot
  # Table1-1: nn.table.q
  # Figure1-2: nn.tsne.q
  # Figure1-3: nn.network.q
  # Figure1-4: nn.ic50.q
  
  
  
  
  ###################################################################################################
  #Module2:Prediction
  #####Upload both files#####
  # example or user data
  dat.source4 <- reactiveValues(source = NULL)
  
  #Example
  observeEvent(input$submitexample.p, {
    showModal(modalDialog(
      size = "m",
      HTML(
        "Analysis started. Please click “Dismiss” to view results. It make take a minute to run."
      ),
      easyClose = TRUE
    ))
    
    dat.source4$source <- "example"
  })
  
  #Message for submit
  observeEvent(input$submitmodel.p, {
    showModal(modalDialog(
      size = "m",
      HTML(
        "Analysis started. Please click “Dismiss” to view results. It make take a minute to run."
      ),
      easyClose = TRUE
    ))
    
    dat.source4$source <- "user"
  })
  
  
  #Read expression data
  observeEvent(input$exp_data.p, {
    exp <- input$exp_data.p
    exp_data <- read.delim(exp$datapath)
    
    if (class(exp_data[1, 1]) != "character") {
      session$sendCustomMessage(type = "resetFileInputHandler", "exp_data.p")
      showModal(modalDialog(
        size = "m",
        HTML(
          "Please ensure that the expression data is the text file with tab-delimited foramt. <br>
        In addition, the format of file is
        the column being the samples and rows being the gene symbols. The first column should be the gene symbols."
        ),
        easyClose = TRUE
      ))
    } else if (length(colnames(exp_data)[-1]) == 1) {
      colnames(exp_data)[1] <- "Gene"
      exp_data$Gene <- toupper(exp_data$Gene)
      pred_data$exp_data = exp_data
      
      showModal(modalDialog(
        size = "m",
        HTML("The format of expression data is correct!"),
        easyClose = TRUE
      ))
      
    } else{
      showModal(modalDialog(
        size = "m",
        HTML(
          "The format of expression data is correct, but we only accept one sample!!!! "
        ),
        easyClose = TRUE
      ))
      
    }
    
    
  })
  
  #Mutation data by gene list
  observeEvent(input$gene_list.p, {
    genelist = toupper(unlist(strsplit(input$gene_list.p, split = "\n")))
    genelist = genelist[which(genelist != "")]
    
    if (length(genelist) != 0) {
      mutation_table  <- data.frame(Gene = genelist, stringsAsFactors = F)
      genelist.update = ShinyDeepDR_CheckGeneSymbol(
        inputData = mutation_table ,
        data.type = "mut",
        alias.table = "data/ccle_exp_and_mut_with_gene_alias.RData"
      )
      
      if (sum(is.na(genelist.update$gene.symbol.update)) == 0) {
        pred_data$genelist <- genelist.update$gene.symbol.update
        showModal(modalDialog(
          size = "m",
          HTML("All gene symbols are valid."),
          easyClose = TRUE
        ))
        
      } else{
        pred_data$genelist <- genelist.update$gene.symbol.update[which(!is.na( genelist.update$gene.symbol.update))]
        showModal(modalDialog(size = "m",
                              HTML(
                                paste(
                                  "Invalid gene symbols : ",
                                  paste(c(genelist.update$gene.symbol[which(is.na(genelist.update$gene.symbol.update))]), collapse = "; "),
                                  "<br>",
                                  "These genes will be excluded from the analysis."
                                )
                              ),
                              easyClose = TRUE))
      }
      
      
    }
  })
  
  #Read mutation data (MAF files)
  observeEvent(input$mut_data.p, {
    mut <- input$mut_data.p
    mut_data <- read.delim(mut$datapath)
    colnames(mut_data) <- toupper(colnames(mut_data))
    
    
    if (sum(colnames(mut_data) %in% toupper(
      c(
        "Hugo_Symbol" ,
        "Variant_Classification",
        "Tumor_Sample_Barcode"
      )
    )) == 3) {
      mut_data = mut_data[, toupper(c(
        "Hugo_Symbol" ,
        "Variant_Classification",
        "Tumor_Sample_Barcode"
      ))]
      mut_data$HUGO_SYMBOL <- toupper(mut_data$HUGO_SYMBOL)
      pred_data$mut_data = mut_data
      session$sendCustomMessage(type = "resetFileInputHandler", "mut_data.p")
      
      if (length(unique(mut_data$TUMOR_SAMPLE_BARCODE)) != 1) {
        showModal(modalDialog(
          size = "m",
          HTML("Multiple samples!! <br> We only accept one sample!!"),
          easyClose = TRUE
          
        ))
        
      } else{
        showModal(modalDialog(
          size = "m",
          HTML(
            "If you already pasted a gene list, you do not need to submit the MAF file. Otherwise, the MAF file will take precedence the list."
          ),
          easyClose = TRUE
        ))
        
        
      }
      
      
    } else{
      session$sendCustomMessage(type = "resetFileInputHandler", "mut_data.p")
      showModal(modalDialog(
        size = "m",
        HTML(
          "Wrong format.<br>
        Please ensure that the mutation data is a MAF file or text file with Hugo_Symbol, Variant_Classification and
             Tumor_Sample_Barcode columns."
        ),
        easyClose = TRUE
      ))
      
      
    }
    
  })
  
  #Prediction both
  
  # Table2-1: pd.table.both
  pd_table_both <- reactive({
    if (is.null(dat.source4$source))
      return()
    
    if (dat.source4$source ==  "example") {
      exp_data =  ex.exp
      mut_data =  ex.mut
      genelist = NA
    }
    
    if (dat.source4$source ==  "user") {
      exp_data =  pred_data$exp_data
      mut_data =  pred_data$mut_data
      
      if(!is.data.frame(mut_data)){
        genelist = pred_data$genelist
        
        print("Use genelist")
      }
      
      if(is.data.frame(mut_data)){
        genelist = NA
        print("Use mutation")
      }
      
    }
    
    data.values <- input$exp_type.p
    log_transform <- input$log_transform.p
    
    #select compaired data type
    cell_type <- input$cell_type1.p
    
    if (!is.data.frame(mut_data) &&  is.na(genelist)) {
      showModal(modalDialog(
        size = "m",
        HTML("Please upload mutation file or paste gene list!"),
        easyClose = TRUE
      ))
      return()
    }
    
    if (is.data.frame(mut_data)) {
      #Check Sample names
      if (sum(colnames(exp_data)[2:ncol(exp_data)] != unique(mut_data$TUMOR_SAMPLE_BARCODE)) >
          0) {
        showModal(modalDialog(
          size = "m",
          HTML(
            "Sample name is not matched between expresstion and mutation data"
          ),
          easyClose = TRUE
          
        ))
      }
    }
    
    #Expression data
    exp_out <- ShinyDeepDR_input_exp(
      inputData = exp_data,
      data.values = data.values,
      selectID = NA,
      log2_transform = log_transform ,
      alias.table = "data/ccle_exp_and_mut_with_gene_alias.RData"
    )
    
    #Mutation data
    if (is.data.frame(mut_data) && is.na(genelist)) {
      mut_out <- ShinyDeepDR_input_mut(
        inputData = mut_data,
        data.format = "MAF",
        alias.table = "data/ccle_exp_and_mut_with_gene_alias.RData"
      )
    }else {
      mut_out <- ShinyDeepDR_input_mut(
        inputData = genelist,
        data.format = "LIST",
        sample_names = exp_out@sample_ID[1],
        alias.table = "data/ccle_exp_and_mut_with_gene_alias.RData"
      )
    }
    
    exp_data.input <- exp_out@exp_final
    rownames(exp_data.input) <- exp_data.input$Gene
    if(class(exp_data.input[,-1]) =="character"){
    exp_data.input[,2] <- as.numeric(exp_data.input[,2])
    }
    #exp_data.input <- data.frame(exp_data.input[,-1])
    exp_data.input.t <-
      matrix(
        t(exp_data.input[, -1]),
        nrow = ncol(exp_data.input) - 1,
        ncol = nrow(exp_data.input),
        dimnames = list(colnames(exp_data.input)[-1],
                        rownames(exp_data.input))
      )
    
    mut_data.input <- mut_out@mutation_table
    if(class(mut_data.input[,-1]) =="character"){
      mut_data.input[,2] <- as.numeric(mut_data.input[,2])
    }
    rownames(mut_data.input) <- mut_data.input$Gene
    mut_data.input.t <-
      matrix(
        t(mut_data.input[, -1]),
        nrow = ncol(mut_data.input) - 1,
        ncol = nrow(mut_data.input),
        dimnames = list(colnames(mut_data.input)[-1],
                        mut_data.input$Gene)
      )
    
    
    dt <- data.table(Pred_Gen(
      mut_data = mut_data.input.t,
      exp_data = exp_data.input.t,
      model =  model
    ))
    colnames(dt) <- drug_names
    dt <- as.data.frame(t(dt))
    colnames(dt) <-
      "Log_predict_IC50" #rownames(exp_data.input.t)
    dt <-
      dt %>% add_column(DrugName = drug_names, .before =  "Log_predict_IC50")
    drug.info.drugname <- drug.info$colname
    drug.info.drugname[256] <- "Salubrinal"
    dt <- dt[match(dt$DrugName, drug.info.drugname), ]
    
    #add sample information
    table.out <- ShinyDeepDR_predict_drug_table(
      inputData = dt,
      cell_type = cell_type,
      ccle_ic50Pred = ccle_ic50_Pred ,
      gdsc_ic50 = drug.ic50,
      tcga_ic50Pred = tcga_ic50Pred,
      gdsc_drugInfo = drug.info[, c(1:12, 15)]
    )#Pred_performance_FullDeepDR
    
    table.out <-
      table.out[order(table.out$Log_predict_IC50, decreasing = F), ]
    
    table.out
    
  })
  
  #Figure 2-1 table
  output$pd.table.both <- renderDataTable({
    table.out =  pd_table_both()
    cell_type <- input$cell_type1.p
    
    if (is.null(table.out))
      return()
    
    if (cell_type == "CCLE") {
      table.out$PubChem <- ShinyDeepDR_CreateLink(table.out$PubChem_CID)
      table.out$PubChem[which(is.na(table.out$PubChem_CID))] <-
        NA
      
      colnames(table.out) <-
        c(
          "Drug name",
          "Predicted IC50 (log uM)",
          "GDSC IC50 min",
          #hide
          "GDSC IC50 max",
          #hide
          "Percentile in GDSC",
          "Pred IC50 min",
          #hide
          "Pred IC50 max",
          #hide
          "Percentile in prediction",
          "colname",
          #hide
          "Synonyms",
          "Action",
          "Clinical stage",
          "Putative target",
          "Targeted pathway",
          "PubChem CID",
          "Canonical SMILES",
          "filename",
          "PubChem link",
          "Prediction performance (Corr)"
        )
      
      table.out$`GDSC range (n=622)` <-
        paste0("[",
               table.out$`GDSC IC50 min`,
               " , ",
               table.out$`GDSC IC50 max`,
               "]")
      table.out$`Prediciton range (n=622)` <-
        paste0("[",
               table.out$`Pred IC50 min`,
               " , ",
               table.out$`Pred IC50 max`,
               "]")
      table.out$`Prediction performance (Corr)` <-
        round(table.out$`Prediction performance (Corr)`, digits = 2)
      table.out.final <- table.out[, c(
        "Drug name",
        "Predicted IC50 (log uM)",
        "Prediction performance (Corr)",
        "GDSC range (n=622)",
        "Percentile in GDSC",
        "Prediciton range (n=622)",
        "Percentile in prediction",
        "Synonyms",
        "Action",
        "Clinical stage",
        "Putative target",
        "Targeted pathway",
        "Canonical SMILES",
        "PubChem CID",
        "PubChem link"
      )]
      
    }
    
    if (cell_type == "TCGA") {
      colnames(table.out) <- c(
        "Drug name",
        "Predicted IC50 (log uM)",
        "TCGA-PredIC50 min",
        "TCGA-PredIC50 max",
        "Percentile in prediction",
        "colname",
        "Synonyms",
        "Action",
        "Clinical stage",
        "Putative target",
        "Targeted pathway",
        "PubChem CID",
        "Canonical SMILES",
        "filename",
        "PubChem link",
        "Prediction performance (Corr)"
      )
      
      table.out$`Prediciton range (n=9050)` <-
        paste0(
          "[",
          round(table.out$`TCGA-PredIC50 min`, digits = 2),
          " , ",
          round(table.out$`TCGA-PredIC50 max`, digits = 2),
          "]"
        )
      table.out$`Prediction performance (Corr)` <-
        round(table.out$`Prediction performance (Corr)`, digits = 2)
      
      table.out.final <- table.out[, c(
        "Drug name",
        "Predicted IC50 (log uM)",
        "Prediction performance (Corr)",
        "Prediciton range (n=9050)",
        "Percentile in prediction",
        "Synonyms",
        "Action",
        "Clinical stage",
        "Putative target",
        "Targeted pathway",
        "Canonical SMILES",
        "PubChem CID",
        "PubChem link"
      )]
    }
    
    #return(table.out.final)
    DT::datatable(
     table.out.final[,-ncol(table.out.final)]
      
    ,
    extensions = c('Scroller', 'FixedColumns'), rownames = FALSE,
    options = list(
      pageLength = 10,
      lengthMenu = c(10, 15, 20),
      dom = 'Brtip',
      scrollX = TRUE
      # extensions = c('Buttons', 'FixedColumns'), rownames = FALSE,width = width.t,height = height.t,
      # options = list(
      #     dom = 'Bfrtip',
      #     buttons = c('excel', 'csv'),
      #     scrollX = TRUE
    ),
    selection = list(mode = "single"))
    
    
  })
  
  # Figure2-2:pd.waterfall.both
  output$pd.waterfall.both <- renderPlotly({
    table.out = pd_table_both()
    
    if (is.null(table.out))
      return()
    
    ShinyDeepDR_WaterFall(inputData = table.out)
    
  })
  
  
  
  # Figure2-3:pd.structure.both
    output$pd.structure.both<- renderImage({
      
      table.out = pd_table_both()
    if (is.null(table.out))
      return()
      
      mypngwidth <- session$clientData$pd.structure.both_width
      mypngheight <- session$clientData$pd.structure.both_height
      
      idx =input$pd.table.both_rows_selected
      
      if(is.null(idx)){
        img.file <- "data/GDSC_paper/DrugStructure/Message_plot.png"
        #img <- readPNG(img.file)
        
      }else{
        
        # Generate the png
        img.file <- paste0("data/GDSC_paper/DrugStructure/",table.out$filename[idx],".png")
        #img <- readPNG(img.file)
        
      }
      
      # Return a list containing the filename
      list(src = img.file,
           width = mypngwidth,
           height = mypngheight ,
           alt = "Please select one drug from the table"
      )
    }, deleteFile = FALSE)

     # Text2-3: pd.drug.info.both
    output$pd.drug.info.both<- renderPlot({

      table.out = pd_table_both()
      if (is.null(table.out))
      return()

      
      # mypngwidth <- session$clientData$pd.structure.both_width
      # mypngheight <- session$clientData$pd.structure.both_height
      
      idx =input$pd.table.both_rows_selected
      if(!is.null(idx)){
        
        par(mar = c(0, 0, 0, 0))
        plot(
        NA,
        xlim = c(0, 2),
        ylim = c(0, 4),
        type = "n",
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = "",
        ann = F,
        bty = 'n'
      )
           legend(
          "topleft",
          cex = 1.8,
          legend = c(
            paste(
              " Drug name: ",
              table.out$DrugName[idx],
              "\n",
              "Predicted IC50 (log uM): ",
              table.out$Log_predict_IC50[idx],
              "\n",
              "Putative target: \n",
              table.out$Putative.Target[idx],
              "\n",
              "Targeted pathway: \n",
              table.out$Targeted.process.pathway[idx],
              "\n",
              "PubChem CID: ",
              table.out$PubChem_CID[idx],
              "\n"
            )
          ),
          bty = "n"
        )

      }else{
        par(mar = c(0, 0, 0, 0))
        plot(x = 0:1,                   # Create empty plot
             y = 0:1,
             ann = F,
             bty = "n",
             type = "n",
             xaxt = "n",
             yaxt = "n")
        text(x = 0.5,                   # Add text to empty plot
             y = 0.5,
             "Please select one drug from the table",
             cex = 1.8)
      }
      
    })
  
  # Figure2-4:pd.density.both
  output$pd.density.both <- renderPlotly({
    table.out = pd_table_both()
    cell_type = input$cell_type1.p
    if (is.null(table.out))
      return()
    
    
    idx <- input$pd.table.both_rows_selected
    
    if (is.null(idx)) {
      idx = 1
    }
    
    if (cell_type == "CCLE") {
      ShinyDeepDR_DensityPlot(
        inputData = table.out,
        cell_type = "CCLE",
        select.drug = table.out$DrugName[idx],
        ccle_ic50Pred = ccle_ic50_Pred,
        gdsc_ic50 = drug.ic50
      )
    } else{
      ShinyDeepDR_DensityPlot(
        inputData = table.out,
        cell_type = "TCGA",
        gdsc_ic50 = drug.ic50,
        select.drug = table.out$DrugName[idx],
        tcga_ic50Pred = tcga_ic50Pred
      )
      
    }
  })
  
  # Figure2-5: pd.box.plot.both
  output$pd.box.plot.both <- renderPlotly({
    table.out = pd_table_both()
    cell_type = input$cell_type1.p
    if (is.null(table.out))
      return()
    
    idx <- input$pd.table.both_rows_selected
    
    if (is.null(idx)) {
      idx = 1
    }
    
    if (cell_type == "CCLE") {
      ShinyDeepDR_BoxPlot(
        inputData = table.out,
        cell_type = "CCLE",
        select.drug = table.out$DrugName[idx],
        ccle_ic50Pred = ccle_ic50_Pred[, table.out$colname],
        gdsc_ic50 = drug.ic50,
        gdsc_cellInfo = cell.info
      )
    } else{
      sorted_indicies = match(table.out$DrugName, tcga_ic50Pred[, 1])
      #print(sorted_indicies)
      ShinyDeepDR_BoxPlot(
        inputData = table.out,
        cell_type = "TCGA",
        select.drug = table.out$DrugName[idx],
        tcga_ic50Pred = tcga_ic50Pred[sorted_indicies, ],
        tcga_clinicalData = tcga_clinicalData
      )
      
    }
    
  })
  
  #Figure 2-6 selection cell or cancer types
  observe({
    cell_type = input$cell_type1.p
    
    if (toupper(cell_type) == "CCLE") {
      cell_type_list <- unique(cell.info$primary_disease)
      cell_type_list <-
        cell_type_list <-
        c(sort(cell_type_list[-which(cell_type_list %in% c("Non-Cancerous", "Unknown"))]), "Non-Cancerous", "Unknown")
      
    } else{
      # Can use character(0) to remove all choices
      cell_type_list <- character(0)
    }
    #
    # if (toupper(cell_type) == "TCGA"){
    #   cell_type_list <- unique(tcga_clinicalData$type)
    #   cell_type_list <- sort(cell_type_list)
    # }
    #
    # Can also set the label and select items
    updateSelectInput(
      session = session,
      inputId = "pd.cancer.select.both",
      label = "Select a Cancer Type",
      choices = cell_type_list,
      selected = "Breast Cancer"
    )
    
    
  })
  
  #Figure2-6 pd.box.plot.both
  output$pd.sub.box.plot.both <- renderPlotly({
    table.out = pd_table_both()
    cell_type = input$cell_type1.p
    select.type = input$pd.cancer.select.both
    idx <- input$pd.table.both_rows_selected
    
    if (is.null(table.out))
      return()
    
    if (is.null(idx)) {
      idx = 1
    }
    
    if (cell_type == "CCLE") {
      shinyDeepDR_CCLESubTypeBoxPlot(
        inputData = table.out,
        select.drug = table.out$DrugName[idx],
        ccle_ic50Pred = ccle_ic50_Pred[, table.out$colname],
        gdsc_ic50 = drug.ic50,
        gdsc_cellInfo = cell.info,
        pancan_type = select.type
      )
    } else{
      ShinyDeepDR_EmptyPlot(title = "Cell Line Only")
      
    }
    
  })
  
  # Figure2-7: pd.network.both
  output$pd.network.both <- renderVisNetwork({
    table.out = pd_table_both()
    cell_type = input$cell_type1.p
    
    if (is.null(table.out))
      return()
    
    idx <- input$pd.table.both_rows_selected
    
    if (is.null(idx)) {
      idx = 1
    }
    
    if (cell_type == "CCLE") {
      ShinyDeepDR_network_for_prediction(
        select.drug = table.out$DrugName[idx],
        cell_type = "CCLE",
        gdsc_ic50 = drug.ic50,
        num_nodes = 10
      )
      
      
    } else{
      ShinyDeepDR_network_for_prediction(
        select.drug = table.out$DrugName[idx],
        cell_type = "TCGA",
        gdsc_ic50 =  drug.ic50,
        tcga_ic50Pred = tcga_ic50Pred,
        tcga_clinicalData = tcga_clinicalData,
        num_nodes = 10
      )
      
    }
    
  })
  
  
  ##############################################################################################################
  #Upload mutation only
  # example or user data
  dat.source5 <- reactiveValues(source = NULL)
  #Example
  observeEvent(input$submitexample2.p, {
    showModal(modalDialog(
      size = "m",
      HTML(
        "Analysis started. Please click “Dismiss” to view results. It make take a minute to run."
      ),
      easyClose = TRUE
    ))
    
    dat.source5$source <- "example"
  })
  
  observeEvent(input$submitmodel2.p, {
    showModal(modalDialog(
      size = "m",
      HTML(
        "Analysis started. Please click “Dismiss” to view results. It make take a minute to run."
      ),
      easyClose = TRUE
    ))
    
    dat.source5$source <- "user"
  })
  
  
  #Mutation data by gene list
  observeEvent(input$gene_list2.p, {
    genelist = toupper(unlist(strsplit(input$gene_list2.p, split = "\n")))
    genelist = genelist[which(genelist != "")]
    
    if (length(genelist) != 0) {
      mutation_table  <- data.frame(Gene = genelist, stringsAsFactors = F)
      genelist.update = ShinyDeepDR_CheckGeneSymbol(
        inputData = mutation_table ,
        data.type = "mut",
        alias.table = "data/ccle_exp_and_mut_with_gene_alias.RData"
      )
      
      if (sum(is.na(genelist.update$gene.symbol.update)) == 0) {
        pred_data$genelist2 <- genelist.update$gene.symbol.update
        showModal(modalDialog(
          size = "m",
          HTML("All gene symbols are valid."),
          easyClose = TRUE
        ))
        
      } else{
        pred_data$genelist2 <- genelist.update$gene.symbol.update[which(!is.na( genelist.update$gene.symbol.update))]
        showModal(modalDialog(size = "m",
                              HTML(
                                paste(
                                  "Invalid gene symbols : ",
                                  paste(c(genelist.update$gene.symbol[which(is.na(genelist.update$gene.symbol.update))]), collapse = "; "),
                                  "<br>",
                                  "These genes will be excluded from the analysis."
                                )
                              ),
                              easyClose = TRUE))
      }
      
      
    }
    
  })
  
  #Read mutation data (MAF files)
  observeEvent(input$mut_data2.p, {
    mut <- input$mut_data2.p
    mut_data <- read.delim(mut$datapath)
    colnames(mut_data) <- toupper(colnames(mut_data))
    
    
    if (sum(colnames(mut_data) %in% toupper(
      c(
        "Hugo_Symbol" ,
        "Variant_Classification",
        "Tumor_Sample_Barcode"
      )
    )) == 3) {
      mut_data = mut_data[, toupper(c(
        "Hugo_Symbol" ,
        "Variant_Classification",
        "Tumor_Sample_Barcode"
      ))]
      mut_data$HUGO_SYMBOL <- toupper(mut_data$HUGO_SYMBOL)
      pred_data$mut_data2 = mut_data
      session$sendCustomMessage(type = "resetFileInputHandler", "mut_data2.p")
      
      if (length(unique(mut_data$TUMOR_SAMPLE_BARCODE)) != 1) {
        showModal(modalDialog(
          size = "m",
          HTML("Multiple samples!! <br> We only accept one sample!!"),
          easyClose = TRUE
          
        ))
        
      } else{
        showModal(modalDialog(
          size = "m",
          HTML(
            "If you already pasted a gene list, you do not need to submit the MAF file. Otherwise, the MAF file will take precedence the list."
          ),
          easyClose = TRUE
        ))
        
        
      }
      
      
    } else{
      session$sendCustomMessage(type = "resetFileInputHandler", "mut_data2.p")
      showModal(modalDialog(
        size = "m",
        HTML(
          "Wrong format.<br>
        Please ensure that the mutation data is a MAF file or text file with Hugo_Symbol, Variant_Classification and
             Tumor_Sample_Barcode columns."
        ),
        easyClose = TRUE
      ))
      
      
    }
    
  })
  
  #Prediction mutation only
  # Table2-1: pd.table.mut
  pd_table_mut <- reactive({
    if (is.null(dat.source5$source))
      return()
    
    if (dat.source5$source ==  "example") {
      mut_data =  ex.mut
      genelist = NA
    }
    
    if (dat.source5$source ==  "user") {
      mut_data =  pred_data$mut_data2
      
      if(!is.data.frame(mut_data)){
        genelist = pred_data$genelist2
        
        print("Use genelist")
      }
      
      if(is.data.frame(mut_data)){
        genelist = NA
        print("Use mutation")
      }
    }
    
    cell_type <- input$cell_type2.p
    
    #Mutation data
    if (is.data.frame(mut_data) && is.na(genelist)) {
      mut_out <- ShinyDeepDR_input_mut(
        inputData = mut_data,
        data.format = "MAF",
        alias.table = "data/ccle_exp_and_mut_with_gene_alias.RData"
      )
    } else{
      mut_out <- ShinyDeepDR_input_mut(
        inputData = genelist,
        data.format = "LIST",
        sample_names = "User_Sample",
        alias.table = "data/ccle_exp_and_mut_with_gene_alias.RData"
      )
    }
    
    mut_data.input <- mut_out@mutation_table
    if(class(mut_data.input[,-1]) =="character"){
      mut_data.input[,2] <- as.numeric(mut_data.input[,2])
    }
    
    rownames(mut_data.input) <- mut_data.input$Gene
    mut_data.input.t <-
      matrix(
        t(mut_data.input[, -1]),
        nrow = ncol(mut_data.input) - 1,
        ncol = nrow(mut_data.input),
        dimnames = list(colnames(mut_data.input)[-1],
                        mut_data.input$Gene)
      )
    
    
    dt <- data.table(Pred_Gen2(mut_data = mut_data.input.t,
                               model =  model3))
    colnames(dt) <- drug_names
    dt <- as.data.frame(t(dt))
    colnames(dt) <-
      "Log_predict_IC50" #rownames(exp_data.input.t)
    dt <-
      dt %>% add_column(DrugName = drug_names, .before =  "Log_predict_IC50")
    drug.info.drugname <- drug.info$colname
    drug.info.drugname[256] <- "Salubrinal"
    dt <- dt[match(dt$DrugName, drug.info.drugname), ]
    
    #add sample information
    table.out <-
      ShinyDeepDR_predict_drug_table(
        inputData = dt,
        cell_type = cell_type,
        ccle_ic50Pred = ccle_ic50_Pred ,
        gdsc_ic50 = drug.ic50,
        tcga_ic50Pred = tcga_ic50Pred,
        gdsc_drugInfo = drug.info[, c(1:12, 13)]
      )#Pred_performance_MutDeepDR
    
    
    
    table.out <-
      table.out[order(table.out$Log_predict_IC50, decreasing = F), ]
    table.out
    
  })
  
  output$pd.table.mut <- renderDataTable({
    table.out =  pd_table_mut()
    cell_type <- input$cell_type2.p
    
    if (is.null(table.out))
      return()
    
    if (cell_type == "CCLE") {
      table.out$PubChem <- ShinyDeepDR_CreateLink(table.out$PubChem_CID)
      table.out$PubChem[which(is.na(table.out$PubChem_CID))] <- NA
      
      colnames(table.out) <-
        c(
          "Drug name",
          "Predicted IC50 (log uM)",
          "GDSC IC50 min",
          #hide
          "GDSC IC50 max",
          #hide
          "Percentile in GDSC",
          "Pred IC50 min",
          #hide
          "Pred IC50 max",
          #hide
          "Percentile in prediction",
          "colname",
          #hide
          "Synonyms",
          "Action",
          "Clinical stage",
          "Putative target",
          "Targeted pathway",
          "PubChem CID",
          "Canonical SMILES",
          "filename",
          "PubChem link",
          "Prediction performance (Corr)"
        )
      
      table.out$`GDSC range (n=622)` <-
        paste0("[",
               table.out$`GDSC IC50 min`,
               " , ",
               table.out$`GDSC IC50 max`,
               "]")
      table.out$`Prediciton range (n=622)` <-
        paste0("[",
               table.out$`Pred IC50 min`,
               " , ",
               table.out$`Pred IC50 max`,
               "]")
      table.out$`Prediction performance (Corr)` <-
        round(table.out$`Prediction performance (Corr)`, digits = 2)
      table.out.final <- table.out[, c(
        "Drug name",
        "Predicted IC50 (log uM)",
        "Prediction performance (Corr)",
        "GDSC range (n=622)",
        "Percentile in GDSC",
        "Prediciton range (n=622)",
        "Percentile in prediction",
        "Synonyms",
        "Action",
        "Clinical stage",
        "Putative target",
        "Targeted pathway",
        "Canonical SMILES",
        "PubChem CID",
        "PubChem link"
      )]
      
      
    }
    
    if (cell_type == "TCGA") {
      colnames(table.out) <- c(
        "Drug name",
        "Predicted IC50 (log uM)",
        "TCGA-PredIC50 min",
        "TCGA-PredIC50 max",
        "Percentile in prediction",
        "colname",
        "Synonyms",
        "Action",
        "Clinical stage",
        "Putative target",
        "Targeted pathway",
        "PubChem CID",
        "Canonical SMILES",
        "filename",
        "PubChem link",
        "Prediction performance (Corr)"
      )
      
      table.out$`Prediciton range (n=9050)` <-
        paste0(
          "[",
          round(table.out$`TCGA-PredIC50 min`, digits = 2),
          " , ",
          round(table.out$`TCGA-PredIC50 max`, digits = 2),
          "]"
        )
      table.out$`Prediction performance (Corr)` <-
        round(table.out$`Prediction performance (Corr)`, digits = 2)
      
      table.out.final <- table.out[, c(
        "Drug name",
        "Predicted IC50 (log uM)",
        "Prediction performance (Corr)",
        "Prediciton range (n=9050)",
        "Percentile in prediction",
        "Synonyms",
        "Action",
        "Clinical stage",
        "Putative target",
        "Targeted pathway",
        "Canonical SMILES",
        "PubChem CID",
        "PubChem link"
      )]
    }
    
    
    DT::datatable(
      table.out.final[,-ncol(table.out.final)]
      
    ,
    extensions = c('Scroller', 'FixedColumns'), rownames = FALSE,
    options = list(
      pageLength = 10,
      lengthMenu = c(10, 15, 20),
      dom = 'Brtip',
      scrollX = TRUE
      # extensions = c('Buttons', 'FixedColumns'), rownames = FALSE,width = width.t,height = height.t,
      # options = list(
      #     dom = 'Bfrtip',
      #     buttons = c('excel', 'csv'),
      #     scrollX = TRUE
    ),
    selection = list(mode = "single"))
    
  }, escape = FALSE)
  
  
  # Figure2-2:pd.waterfall.mut
  output$pd.waterfall.mut <- renderPlotly({
    table.out = pd_table_mut()
    if (is.null(table.out))
      return()
    
    ShinyDeepDR_WaterFall(inputData = table.out)
  })
  
 
  # Figure2-3:pd.structure.mut
    output$pd.structure.mut<- renderImage({
      table.out = pd_table_mut()
    
    if (is.null(table.out))
      return()
      
      mypngwidth <- session$clientData$pd.structure.mut_width
      mypngheight <- session$clientData$pd.structure.mut_height
      
      idx =input$pd.table.mut_rows_selected
      
      if(is.null(idx)){
        img.file <- "data/GDSC_paper/DrugStructure/Message_plot.png"
        #img <- readPNG(img.file)
        
      }else{
        
        # Generate the png
        img.file <- paste0("data/GDSC_paper/DrugStructure/",table.out$filename[idx],".png")
        #img <- readPNG(img.file)
        
      }
      
      # Return a list containing the filename
      list(src = img.file,
           width = mypngwidth,
           height = mypngheight ,
           alt = "Please select one drug from the table"
      )
    }, deleteFile = FALSE)

    # Text2-3: pd.drug.info.mut
    output$pd.drug.info.mut<- renderPlot({
      table.out = pd_table_mut()
    
    if (is.null(table.out))
      return()
      
      # mypngwidth <- session$clientData$pd.structure.both_width
      # mypngheight <- session$clientData$pd.structure.both_height
      
      idx =input$pd.table.mut_rows_selected
      if(!is.null(idx)){
        
        par(mar = c(0, 0, 0, 0))
        plot(
        NA,
        xlim = c(0, 2),
        ylim = c(0, 4),
        type = "n",
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = "",
        ann = F,
        bty = 'n'
      )
      legend(
          "topleft",
          cex = 1.8,
          legend = c(
            paste(
              " Drug name: ",
              table.out$DrugName[idx],
              "\n",
              "Predicted IC50 (log uM): ",
              table.out$Log_predict_IC50[idx],
              "\n",
              "Putative target: \n",
              table.out$Putative.Target[idx],
              "\n",
              "Targeted pathway: \n",
              table.out$Targeted.process.pathway[idx],
              "\n",
              "PubChem CID: ",
              table.out$PubChem_CID[idx],
              "\n"
            )
          ),
          bty = "n"
        )

      }else{
        par(mar = c(0, 0, 0, 0))
        plot(x = 0:1,                   # Create empty plot
             y = 0:1,
             ann = F,
             bty = "n",
             type = "n",
             xaxt = "n",
             yaxt = "n")
        text(x = 0.5,                   # Add text to empty plot
             y = 0.5,
             "Please select one drug from the table",
             cex = 1.8)
      }
      
    })
  
  # Figure2-4:pd.density.mut
  output$pd.density.mut <- renderPlotly({
    table.out = pd_table_mut()
    cell_type = input$cell_type2.p
    
    if (is.null(table.out))
      return()
    
    idx <- input$pd.table.mut_rows_selected
    
    if (is.null(idx)) {
      idx = 1
    }
    
    if (cell_type == "CCLE") {
      ShinyDeepDR_DensityPlot(
        inputData = table.out,
        cell_type = "CCLE",
        select.drug = table.out$DrugName[idx],
        ccle_ic50Pred = ccle_ic50_Pred,
        gdsc_ic50 = drug.ic50
      )
    } else{
      ShinyDeepDR_DensityPlot(
        inputData = table.out,
        cell_type = "TCGA",
        gdsc_ic50 = drug.ic50,
        select.drug = table.out$DrugName[idx],
        tcga_ic50Pred = tcga_ic50Pred
      )
      
    }
  })
  
  # Figure2-5: pd.box.plot.mut
  output$pd.box.plot.mut <- renderPlotly({
    table.out = pd_table_mut()
    cell_type = input$cell_type2.p
    
    if (is.null(table.out))
      return()
    
    idx <- input$pd.table.mut_rows_selected
    
    if (is.null(idx)) {
      idx = 1
    }
    
    if (cell_type == "CCLE") {
      ShinyDeepDR_BoxPlot(
        inputData = table.out,
        cell_type = "CCLE",
        select.drug = table.out$DrugName[idx],
        ccle_ic50Pred = ccle_ic50_Pred[, table.out$colname],
        gdsc_ic50 = drug.ic50,
        gdsc_cellInfo = cell.info
      )
    } else{
      sorted_indicies = match(table.out$DrugName, tcga_ic50Pred[, 1])
      #print(sorted_indicies)
      ShinyDeepDR_BoxPlot(
        inputData = table.out,
        cell_type = "TCGA",
        select.drug = table.out$DrugName[idx],
        tcga_ic50Pred = tcga_ic50Pred[sorted_indicies, ],
        tcga_clinicalData = tcga_clinicalData
      )
      
    }
  })
  
  #Figure 2-6 selection cell or cancer types
  observe({
    cell_type = input$cell_type2.p
    
    if (toupper(cell_type) == "CCLE") {
      cell_type_list <- unique(cell.info$primary_disease)
      cell_type_list <-
        cell_type_list <-
        c(sort(cell_type_list[-which(cell_type_list %in% c("Non-Cancerous", "Unknown"))]), "Non-Cancerous", "Unknown")
      
    } else{
      # Can use character(0) to remove all choices
      cell_type_list <- character(0)
    }
    #
    # if (toupper(cell_type) == "TCGA"){
    #   cell_type_list <- unique(tcga_clinicalData$type)
    #   cell_type_list <- sort(cell_type_list)
    # }
    #
    # Can also set the label and select items
    updateSelectInput(
      session = session,
      inputId = "pd.cancer.select.mut",
      label = "Select a Cancer Type",
      choices = cell_type_list,
      selected = "Breast Cancer"
    )
    
  })
  
  #Figure2-6 pd.box.plot.mut
  output$pd.sub.box.plot.mut <- renderPlotly({
    table.out = pd_table_mut()
    cell_type = input$cell_type2.p
    select.type = input$pd.cancer.select.mut
    idx <- input$pd.table.mut_rows_selected
    
    if (is.null(table.out))
      return()
    
    if (is.null(idx)) {
      idx = 1
    }
    
    if (cell_type == "CCLE") {
      shinyDeepDR_CCLESubTypeBoxPlot(
        inputData = table.out,
        select.drug = table.out$DrugName[idx],
        ccle_ic50Pred = ccle_ic50_Pred[, table.out$colname],
        gdsc_ic50 = drug.ic50,
        gdsc_cellInfo = cell.info,
        pancan_type = select.type
      )
    } else{
      ShinyDeepDR_EmptyPlot(title = "Cell Line Only")
      
    }
  })
  
  # Figure2-7: pd.network.mut
  output$pd.network.mut <- renderVisNetwork({
    table.out = pd_table_mut()
    cell_type = input$cell_type2.p
    
    if (is.null(table.out))
      return()
    
    idx <- input$pd.table.mut_rows_selected
    
    if (is.null(idx)) {
      idx = 1
    }
    
    if (cell_type == "CCLE") {
      ShinyDeepDR_network_for_prediction(
        select.drug = table.out$DrugName[idx],
        cell_type = cell_type,
        gdsc_ic50 = drug.ic50,
        num_nodes = 10
      )
      
      
    } else{
      ShinyDeepDR_network_for_prediction(
        select.drug = table.out$DrugName[idx],
        cell_type = "TCGA",
        gdsc_ic50 = drug.ic50,
        tcga_ic50Pred = tcga_ic50Pred,
        tcga_clinicalData = tcga_clinicalData,
        num_nodes = 10
      )
      
    }
  })
  
  ##############################################################################################################
  # example or user data
  dat.source6 <- reactiveValues(source = NULL)
  
  #Example
  observeEvent(input$submitexample3.p, {
    showModal(modalDialog(
      size = "m",
      HTML(
        "Analysis started. Please click “Dismiss” to view results. It make take a minute to run."
      ),
      easyClose = TRUE
    ))
    
    dat.source6$source <- "example"
    
  })
  
  #Upload exp data only
  observeEvent(input$submitmodel3.p, {
    showModal(modalDialog(
      size = "m",
      HTML(
        "Analysis started. Please click “Dismiss” to view results. It make take a minute to run."
      ),
      easyClose = TRUE
    ))
    dat.source6$source <- "user"
  })
  
  #Read expression data
  observeEvent(input$exp_data2.p, {
    exp <- input$exp_data2.p
    exp_data <- read.delim(exp$datapath)
    
    if (class(exp_data[1, 1]) != "character") {
      session$sendCustomMessage(type = "resetFileInputHandler", "exp_data2.p")
      showModal(modalDialog(
        size = "m",
        HTML(
          "Please ensure that the expression data is the text file with tab-delimited foramt. <br>
        In addition, the format of file is
        the column being the samples and rows being the gene symbols. The first column should be the gene symbols."
        ),
        easyClose = TRUE
      ))
    } else{
      colnames(exp_data)[1] <- "Gene"
      exp_data$Gene <- toupper(exp_data$Gene)
      pred_data$exp_data2 = exp_data
      showModal(modalDialog(
        size = "m",
        HTML("The format of expression data is correct!"),
        easyClose = TRUE
      ))
    }
    
    
  })
  
  
  #Prediction --> gene expression only
  # Table2-1: pd.table.exp
  pd_table_exp <- reactive({
    if (is.null(dat.source6$source))
      return()
    
    if (dat.source6$source ==  "example") {
      exp_data =  ex.exp
    }
    
    if (dat.source6$source ==  "user") {
      exp_data =  pred_data$exp_data2
    }
    
    
    data.values <- input$exp_type2.p
    log_transform <- input$log_transform2.p
    cell_type <- input$cell_type3.p
    
    
    #Expression data
    exp_out <- ShinyDeepDR_input_exp(
      inputData = exp_data,
      data.values = data.values,
      selectID = NA,
      log2_transform = log_transform ,
      alias.table = "data/ccle_exp_and_mut_with_gene_alias.RData"
    )
    
    exp_data.input <- exp_out@exp_final
    if(class(exp_data.input[,-1]) =="character"){
      exp_data.input[,2] <- as.numeric(exp_data.input[,2])
    }
    rownames(exp_data.input) <- exp_data.input$Gene
    #exp_data.input <- data.frame(exp_data.input[,-1])
    exp_data.input.t <-
      matrix(
        t(exp_data.input[, -1]),
        nrow = ncol(exp_data.input) - 1,
        ncol = nrow(exp_data.input),
        dimnames = list(colnames(exp_data.input)[-1],
                        rownames(exp_data.input))
      )
    
    
    
    dt <- data.table(Pred_Gen3(exp_data = exp_data.input.t,
                               model =  model2))
    colnames(dt) <- drug_names
    dt <- as.data.frame(t(dt))
    colnames(dt) <-  "Log_predict_IC50" #rownames(exp_data.input.t)
    dt <-
      dt %>% add_column(DrugName = drug_names, .before =  "Log_predict_IC50")
    drug.info.drugname <- drug.info$colname
    drug.info.drugname[256] <- "Salubrinal"
    dt <- dt[match(dt$DrugName, drug.info.drugname), ]
    
    #add sample information
    table.out <-
      ShinyDeepDR_predict_drug_table(
        inputData = dt,
        cell_type = cell_type,
        ccle_ic50Pred = ccle_ic50_Pred ,
        gdsc_ic50 = drug.ic50,
        tcga_ic50Pred = tcga_ic50Pred,
        gdsc_drugInfo = drug.info[, c(1:12, 14)]
      )
    
    
    table.out <-
      table.out[order(table.out$Log_predict_IC50, decreasing = F), ]
    table.out
  })
  
  output$pd.table.exp <- renderDataTable({
    table.out =  pd_table_exp()
    cell_type <- input$cell_type3.p
    
    if (is.null(table.out))
      return()
    
    if (cell_type == "CCLE") {
      table.out$PubChem <- ShinyDeepDR_CreateLink(table.out$PubChem_CID)
      table.out$PubChem[which(is.na(table.out$PubChem_CID))] <- NA
      
      colnames(table.out) <-
        c(
          "Drug name",
          "Predicted IC50 (log uM)",
          "GDSC IC50 min",
          #hide
          "GDSC IC50 max",
          #hide
          "Percentile in GDSC",
          "Pred IC50 min",
          #hide
          "Pred IC50 max",
          #hide
          "Percentile in prediction",
          "colname",
          #hide
          "Synonyms",
          "Action",
          "Clinical stage",
          "Putative target",
          "Targeted pathway",
          "PubChem CID",
          "Canonical SMILES",
          "filename",
          "PubChem link",
          "Prediction performance (Corr)"
        )
      
      table.out$`GDSC range (n=622)` <-
        paste0("[",
               table.out$`GDSC IC50 min`,
               " , ",
               table.out$`GDSC IC50 max`,
               "]")
      table.out$`Prediciton range (n=622)` <-
        paste0("[",
               table.out$`Pred IC50 min`,
               " , ",
               table.out$`Pred IC50 max`,
               "]")
      table.out$`Prediction performance (Corr)` <-
        round(table.out$`Prediction performance (Corr)`, digits = 2)
      table.out.final <- table.out[, c(
        "Drug name",
        "Predicted IC50 (log uM)",
        "Prediction performance (Corr)",
        "GDSC range (n=622)",
        "Percentile in GDSC",
        "Prediciton range (n=622)",
        "Percentile in prediction",
        "Synonyms",
        "Action",
        "Clinical stage",
        "Putative target",
        "Targeted pathway",
        "Canonical SMILES",
        "PubChem CID",
        "PubChem link"
      )]
      
      
    }
    
    if (cell_type == "TCGA") {
      colnames(table.out) <- c(
        "Drug name",
        "Predicted IC50 (log uM)",
        "TCGA-PredIC50 min",
        "TCGA-PredIC50 max",
        "Percentile in prediction",
        "colname",
        "Synonyms",
        "Action",
        "Clinical stage",
        "Putative target",
        "Targeted pathway",
        "PubChem CID",
        "Canonical SMILES",
        "filename",
        "PubChem link",
        "Prediction performance (Corr)"
      )
      
      table.out$`Prediciton range (n=9050)` <-
        paste0(
          "[",
          round(table.out$`TCGA-PredIC50 min`, digits = 2),
          " , ",
          round(table.out$`TCGA-PredIC50 max`, digits = 2),
          "]"
        )
      table.out$`Prediction performance (Corr)` <-
        round(table.out$`Prediction performance (Corr)`, digits = 2)
      
      table.out.final <- table.out[, c(
        "Drug name",
        "Predicted IC50 (log uM)",
        "Prediction performance (Corr)",
        "Prediciton range (n=9050)",
        "Percentile in prediction",
        "Synonyms",
        "Action",
        "Clinical stage",
        "Putative target",
        "Targeted pathway",
        "Canonical SMILES",
        "PubChem CID",
        "PubChem link"
      )]
      
    }
    
    DT::datatable(
      table.out.final[,-ncol(table.out.final)]
    ,
    extensions = c('Scroller', 'FixedColumns'),
    rownames = FALSE,
    options = list(
      pageLength = 10,
      lengthMenu = c(10, 15, 20),
      dom = 'Brtip',
      scrollX = TRUE
      # extensions = c('Buttons', 'FixedColumns'),
      # options = list(
      #     dom = 'Bfrtip',
      #     buttons = c('excel', 'csv'),
      #     scrollX = TRUE
    ),
    selection = list(mode = "single"))
    
  }, escape = TRUE)
  
  
  # Figure2-2:pd.waterfall.exp
  output$pd.waterfall.exp <- renderPlotly({
    table.out = pd_table_exp()
    if (is.null(table.out))
      return()
    
    ShinyDeepDR_WaterFall(inputData = table.out)
  })
  
# Figure2-3:pd.structure.exp
    output$pd.structure.exp<- renderImage({
    table.out = pd_table_exp()
    
    if (is.null(table.out))
      return()
      
      mypngwidth <- session$clientData$pd.structure.exp_width
      mypngheight <- session$clientData$pd.structure.exp_height
      
      idx =input$pd.table.exp_rows_selected
      
      if(is.null(idx)){
        img.file <- "data/GDSC_paper/DrugStructure/Message_plot.png"
        #img <- readPNG(img.file)
        
      }else{
        
        # Generate the png
        img.file <- paste0("data/GDSC_paper/DrugStructure/",table.out$filename[idx],".png")
        #img <- readPNG(img.file)
        
      }
      
      # Return a list containing the filename
      list(src = img.file,
           width = mypngwidth,
           height = mypngheight ,
           alt = "Please select one drug from the table"
      )
    }, deleteFile = FALSE)

    # Text2-3: pd.drug.info.exp
    output$pd.drug.info.exp<- renderPlot({
      table.out = pd_table_exp()
    
    if (is.null(table.out))
      return()
      
      idx =input$pd.table.exp_rows_selected
      if(!is.null(idx)){
        
        
        par(mar = c(0, 0, 0, 0))
        plot(
        NA,
        xlim = c(0, 2),
        ylim = c(0, 4),
        type = "n",
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = "",
        ann = F,
        bty = 'n'
      )
           legend(
          "topleft",
          cex = 1.8,
          legend = c(
            paste(
              " Drug name: ",
              table.out$DrugName[idx],
              "\n",
              "Predicted IC50 (log uM): ",
              table.out$Log_predict_IC50[idx],
              "\n",
              "Putative target: \n",
              table.out$Putative.Target[idx],
              "\n",
              "Targeted pathway: \n",
              table.out$Targeted.process.pathway[idx],
              "\n",
              "PubChem CID: ",
              table.out$PubChem_CID[idx],
              "\n"
            )
          ),
          bty = "n"
        )

      }else{
        par(mar = c(0, 0, 0, 0))
        plot(x = 0:1,                   # Create empty plot
             y = 0:1,
             ann = F,
             bty = "n",
             type = "n",
             xaxt = "n",
             yaxt = "n")
        text(x = 0.5,                   # Add text to empty plot
             y = 0.5,
             "Please select one drug from the table",
             cex = 1.8)
      }
        
      
      
    })
  
  # Figure2-4:pd.density.exp
  output$pd.density.exp <- renderPlotly({
    table.out = pd_table_exp()
    cell_type = input$cell_type3.p
    
    if (is.null(table.out))
      return()
    
    
    idx <- input$pd.table.exp_rows_selected
    
    if (is.null(idx)) {
      idx = 1
    }
    
    if (cell_type == "CCLE") {
      ShinyDeepDR_DensityPlot(
        inputData = table.out,
        cell_type = "CCLE",
        select.drug = table.out$DrugName[idx],
        ccle_ic50Pred = ccle_ic50_Pred,
        gdsc_ic50 = drug.ic50
      )
    } else{
      ShinyDeepDR_DensityPlot(
        inputData = table.out,
        cell_type = "TCGA",
        gdsc_ic50 = drug.ic50,
        select.drug = table.out$DrugName[idx],
        tcga_ic50Pred = tcga_ic50Pred
      )
      
    }
  })
  
  
  # Figure2-5: pd.box.plot.exp
  output$pd.box.plot.exp <- renderPlotly({
    table.out = pd_table_exp()
    cell_type = input$cell_type3.p
    if (is.null(table.out))
      return()
    
    idx <- input$pd.table.exp_rows_selected
    
    if (is.null(idx)) {
      idx = 1
    }
    
    if (cell_type == "CCLE") {
      ShinyDeepDR_BoxPlot(
        inputData = table.out,
        cell_type = "CCLE",
        select.drug = table.out$DrugName[idx],
        ccle_ic50Pred = ccle_ic50_Pred[, table.out$colname],
        gdsc_ic50 = drug.ic50,
        gdsc_cellInfo = cell.info
      )
    } else{
      sorted_indicies = match(table.out$DrugName, tcga_ic50Pred[, 1])
      #print(sorted_indicies)
      ShinyDeepDR_BoxPlot(
        inputData = table.out,
        cell_type = "TCGA",
        select.drug = table.out$DrugName[idx],
        tcga_ic50Pred = tcga_ic50Pred[sorted_indicies, ],
        tcga_clinicalData = tcga_clinicalData
      )
      
    }
  })
  
  #Figure 2-6 selection cell or cancer types
  observe({
    cell_type = input$cell_type3.p
    
    if (toupper(cell_type) == "CCLE") {
      cell_type_list <- unique(cell.info$primary_disease)
      cell_type_list <-
        c(sort(cell_type_list[-which(cell_type_list %in% c("Non-Cancerous", "Unknown"))]), "Non-Cancerous", "Unknown")
      
    } else{
      # Can use character(0) to remove all choices
      cell_type_list <- character(0)
    }
    #
    # if (toupper(cell_type) == "TCGA"){
    #   cell_type_list <- unique(tcga_clinicalData$type)
    #   cell_type_list <- sort(cell_type_list)
    # }
    #
    # Can also set the label and select items
    updateSelectInput(
      session = session,
      inputId = "pd.cancer.select.exp",
      label = "Select a Cancer Type",
      choices = cell_type_list,
      selected = "Breast Cancer"
    )
    
    
    
  })
  
  
  #Figure2-6 pd.box.plot.exp
  output$pd.sub.box.plot.exp <- renderPlotly({
    table.out = pd_table_exp()
    cell_type = input$cell_type3.p
    select.type = input$pd.cancer.select.exp
    idx <- input$pd.table.exp_rows_selected
    
    if (is.null(table.out))
      return()
    
    if (is.null(idx)) {
      idx = 1
    }
    
    if (cell_type == "CCLE") {
      shinyDeepDR_CCLESubTypeBoxPlot(
        inputData = table.out,
        select.drug = table.out$DrugName[idx],
        ccle_ic50Pred = ccle_ic50_Pred[, table.out$colname],
        gdsc_ic50 = drug.ic50,
        gdsc_cellInfo = cell.info,
        pancan_type = select.type
      )
    } else{
      ShinyDeepDR_EmptyPlot(title = "Cell Line Only")
      
    }
    
  })
  
  # Figure2-7: pd.network.exp
  output$pd.network.exp <- renderVisNetwork({
    table.out = pd_table_exp()
    cell_type = input$cell_type3.p
    if (is.null(table.out))
      return()
    
    idx <- input$pd.table.exp_rows_selected
    
    if (is.null(idx)) {
      idx = 1
    }
    
    if (cell_type == "CCLE") {
      ShinyDeepDR_network_for_prediction(
        select.drug = table.out$DrugName[idx],
        cell_type = "CCLE",
        gdsc_ic50 = drug.ic50,
        num_nodes = 10
      )
      
      
    } else{
      
      ShinyDeepDR_network_for_prediction(
        select.drug = table.out$DrugName[idx],
        cell_type = "TCGA",
        gdsc_ic50 =  drug.ic50,
        tcga_ic50Pred = tcga_ic50Pred,
        tcga_clinicalData = tcga_clinicalData,
        num_nodes = 10
      )
      
    }
  })
  
  
  
  #Show description for each plot
  # Table2-1: pd.table.q
  # Figure2-2:pd.waterfall.q
  # Figure2-3:pd.structure.q
  # Text2-3: pd.drug.info.q
  # Figure2-4:pd.density.q
  # Figure2-5: pd.box.plot.q
  # Figure2-6: pd.network.q
  
})
