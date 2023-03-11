#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#' ---
#' title: DepLink: A missing link between cancer gene dependency and drug response
#' author: Y-C Chiu, T. Nayak, L-J Wang, and Y. Chen
#' date: Novemeber 2021
#' Affiliation: Greehey Children's Cancer Research Institute,
#'              UT Health San Antonio, San Antonio, TX
#' ---

#' The shiny server main program including ui() and server()

# Nov. 2021, app.R created, L-J. Wang
# Nov. 19, 2021, Header added.  Y. Chen

# You may need to install following packages
#install.packages("rintrojs")
#install.packages("shinyBS")
#install.packages("visNetwork")
library("htmltools")
library("htmlwidgets")
library("rintrojs")
library("shiny")
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
library("visNetwork")
library("stringr")
library("dplyr")
#library("sqldf")
library("shinycssloaders")
library("shinyjs")
library("igraph")




#TODO: separate the 3 edge components into g2d, d2d, and g2g and make 3 dropdowns. figure out how to do clusters with the edges in R
#TODO: change the size of the genes, add general info about the network in a table. when selecting gene without cluster, just select gene


#setwd("~/Desktop/DepLink_Michael/")


ui <- htmlTemplate(
  "index.html",
  ########
  #Privacy policy:
  #UI ID: privacytabg1
  #If you need multiple privacy ids, please add number, such as “privacytag1”.
  privacytag1 = actionLink(
    inputId = "privacytag1",
    label = "Privacy Policy",
    style = "color:white;border-bottom: 1px solid #a9529e;"
  ),
  
  #Show description for each plot
  # Table1-1: nn.table.q
  # Figure1-3: nn.network.q
  nnboth1 = circleButton(
    inputId = "nnboth1",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  nnboth3 = circleButton(
    inputId = "nnboth3",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  
  nnmut1 = circleButton(
    inputId = "nnmut1",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  nnmut3 = circleButton(
    inputId = "nnmut3",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  
  nnexp1 = circleButton(
    inputId = "nnexp1",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  nnexp3 = circleButton(
    inputId = "nnexp3",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  helptab3p1 = actionLink(inputId = "helptab3p1", label = "Input Data Format"),
  helptab1p1 = actionLink(inputId = "helptab1p1", label = "Find Similar Samples (Find Sample)"),
  helptab2p1 = actionLink(inputId = "helptab2p1", label = "Predict Drug Response (Find Drug)"),
  
  helptab3p2 = actionLink(inputId = "helptab3p2", label = "Input Data Format"),
  helptab1p2 = actionLink(inputId = "helptab1p2", label = "Find Similar Samples (Find Sample)"),
  helptab2p2 = actionLink(inputId = "helptab2p2", label = "Predict Drug Response (Find Drug)"),
  
  helptab3p3 = actionLink(inputId = "helptab3p3", label = "Input Data Format"),
  helptab1p3 = actionLink(inputId = "helptab1p3", label = "Find Similar Samples (Find Sample)"),
  helptab2p3 = actionLink(inputId = "helptab2p3", label = "Predict Drug Response (Find Drug)"),
  
  helptab3p4 = actionLink(inputId = "helptab3p4", label = "Input Data Format"),
  helptab1p4 = actionLink(inputId = "helptab1p4", label = "Find Similar Samples (Find Sample)"),
  helptab2p4 = actionLink(inputId = "helptab2p4", label = "Predict Drug Response (Find Drug)"),
  
  helptab3p5 = actionLink(inputId = "helptab3p5", label = "Input Data Format"),
  helptab1p5 = actionLink(inputId = "helptab1p5", label = "Find Similar Samples (Find Sample)"),
  helptab2p5 = actionLink(inputId = "helptab2p5", label = "Predict Drug Response (Find Drug)"),
  
  help_input1 = actionLink(
    inputId = "help_input1",
    label = "Learn more about input data format and download examples",
    style = "color:blueviolet;
                                         border-bottom: 1px solid #a9529e;"
  ),
  help_input2 = actionLink(
    inputId = "help_input2",
    label = "Learn more about input data format and download examples",
    style = "color:blueviolet;
                                         border-bottom: 1px solid #a9529e;"
  ),
  
  pdboth1 = circleButton(
    inputId = "pdboth1",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  pdboth2 = circleButton(
    inputId = "pdboth2",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  pdboth3 = circleButton(
    inputId = "pdboth3",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  pdboth5 = circleButton(
    inputId = "pdboth5",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  pdboth6 = circleButton(
    inputId = "pdboth6",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  
  #NEW MODAL
  pdboth6sub = circleButton(
    inputId = "pdboth6sub",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  pdboth7 = circleButton(
    inputId = "pdboth7",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  
  pdmut1 = circleButton(
    inputId = "pdmut1",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  pdmut2 = circleButton(
    inputId = "pdmut2",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  pdmut3 = circleButton(
    inputId = "pdmut3",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  pdmut5 = circleButton(
    inputId = "pdmut5",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  pdmut6 = circleButton(
    inputId = "pdmut6",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  #NEW MODAL
  pdmut6sub = circleButton(
    inputId = "pdmut6sub",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  pdmut7 = circleButton(
    inputId = "pdmut7",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  
  pdexp1 = circleButton(
    inputId = "pdexp1",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  pdexp2 = circleButton(
    inputId = "pdexp2",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  pdexp3 = circleButton(
    inputId = "pdexp3",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  pdexp5 = circleButton(
    inputId = "pdexp5",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  pdexp6 = circleButton(
    inputId = "pdexp6",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  #NEW MODAL
  pdexp6sub = circleButton(
    inputId = "pdexp6sub",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  pdexp7 = circleButton(
    inputId = "pdexp7",
    icon = icon("question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  
  #Module1:Nearest Neighbors:
  #####Upload both files#####
  submitmodel = actionBttn(
    inputId = "submitmodel",
    label = "Submit",
    size = "md",
    icon = icon("play-circle"),
    style = "material-flat",
    color = "royal"
  ),
  
  submitexample = actionBttn(
    inputId = "submitexample",
    label = "Run Example",
    size = "md",
    icon = icon("rocket"),
    style = "material-flat",
    color = "primary"
  ),
  #example ID--> example.nn.both
  
  mut_upload = fileInput(
    "mut_data",
    NULL,
    multiple = FALSE,
    accept = c("text/csv",".maf",
               "text/comma-separated-values,text/plain",
               ".csv"),
    width = "100%"
  ),
  
  gene_list = textAreaInput(
    "gene_list",
    label = "",
    width = "100%",
    height = "200px",
    placeholder = "Paste gene symbols here as separate lines, for example:\nCDK4\nCDK6\nCCND1"
  ),
  exp_upload = fileInput(
    "exp_data",
    NULL,
    multiple = FALSE,
    accept = c("text/csv",
               "text/comma-separated-values,text/plain",
               ".csv"),
    width = "100%"
  ),
  
  exp_type = prettyRadioButtons(
    "exp_type",
    label = "Data Type",
    choices = c("FPKM", "TPM"),
    inline = TRUE,
    shape = "square",
    selected =  "TPM"
  ),
  
  log_transform = prettyRadioButtons(
    "log_transform",
    label = "Log Transformed (log2(data+1))?",
    choiceNames = c("Yes", "No"),
    choiceValues = c(1, 2),
    inline = TRUE,
    shape = "square",
    selected = 1
  ),
  
  
  cell_type = prettyRadioButtons(
    inputId = "cell_type1",
    label = "",
    choiceNames = c("Cell line", "Tumor"),
    choiceValues = c("CCLE", "TCGA"),
    inline = TRUE,
    shape = "square",
    selected = "CCLE"
  ),
  
  #Show results:
  # Table1-1: nn.table.both
  # Figure1-3: nn.network.both
  
  nn.table.both = shinycssloaders::withSpinner(dataTableOutput("nn.table.both")),
  nn.network.both = shinycssloaders::withSpinner(visNetworkOutput("nn.network.both")),
  
  
  
  #####Upload mutation only#####
  submitmodel2 = actionBttn(
    inputId = "submitmodel2",
    label = "Submit",
    size = "md",
    icon = icon("play-circle"),
    style = "material-flat",
    color = "royal"
  ),
  submitexample2 = actionBttn(
    inputId = "submitexample2",
    label = "Run Example",
    size = "md",
    icon = icon("rocket"),
    style = "material-flat",
    color = "primary"
  ),
  #example ID -->example.nn.mut
  
  mut_upload2 = fileInput(
    "mut_data2",
    NULL,
    multiple = FALSE,
    accept = c("text/csv",".maf",
               "text/comma-separated-values,text/plain",
               ".csv"),
    width = "100%"
  ),
  
  gene_list2 = textAreaInput(
    "gene_list2",
    label = "",
    width = "100%",
    height = "200px",
    placeholder = "Paste gene symbols here as separate lines, for example:\nCDK4\nCDK6\nCCND1"
  ),
  cell_type2 = prettyRadioButtons(
    "cell_type2",
    label = "",
    choiceNames = c("Cell line", "Tumor"),
    choiceValues = c("CCLE", "TCGA"),
    inline = TRUE,
    shape = "square",
    selected = "CCLE"
  ),
  #Show results:
  # Table1-1: nn.table.mut
  # Figure1-3: nn.network.mut
  
  nn.table.mut = shinycssloaders::withSpinner(dataTableOutput("nn.table.mut")),
  nn.network.mut = shinycssloaders::withSpinner(visNetworkOutput("nn.network.mut")),
  
  
  #####upload gene expression only#####
  submitmodel3 = actionBttn(
    inputId = "submitmodel3",
    label = "Submit",
    size = "md",
    icon = icon("play-circle"),
    style = "material-flat",
    color = "royal"
  ),
  submitexample3 = actionBttn(
    inputId = "submitexample3",
    label = "Run Example",
    size = "md",
    icon = icon("rocket"),
    style = "material-flat",
    color = "primary"
  ),
  
  #example ID -->example.nn.exp
  
  exp_upload2 = fileInput(
    "exp_data2",
    NULL,
    multiple = FALSE,
    accept = c("text/csv",
               "text/comma-separated-values,text/plain",
               ".csv"),
    width = "100%"
  ),
  exp_type2 = prettyRadioButtons(
    "exp_type2",
    label = "Data Type",
    choices =
      c("FPKM", "TPM"),
    inline = TRUE,
    shape = "square",
    selected = "TPM"
  ),
  log_transform2 = prettyRadioButtons(
    "log_transform2",
    label = "Log Transformed (log2(data+1))?",
    choiceNames = c("Yes", "No"),
    choiceValues = c(1, 2),
    inline = TRUE,
    shape = "square",
    selected = 1
  ),
  cell_type3 = prettyRadioButtons(
    "cell_type3",
    label = "",
    choiceNames = c("Cell line",  "Tumor"),
    choiceValues = c("CCLE", "TCGA"),
    inline = TRUE,
    shape = "square",
    selected = "CCLE"
  ),
  #Show results:
  # Table1-1: nn.table.exp
  # Figure1-3: nn.network.exp
  
  
  nn.table.exp = shinycssloaders::withSpinner(dataTableOutput("nn.table.exp")),
  nn.network.exp = shinycssloaders::withSpinner(visNetworkOutput("nn.network.exp")),
  
  #Show description for each plot
  # Table1-1: nn.table.q
  # Figure1-3: nn.network.q
  # nnboth1 = circleButton(
  #   inputId = "nnboth1",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # nnboth3 = circleButton(
  #   inputId = "nnboth3",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # 
  # nnmut1 = circleButton(
  #   inputId = "nnmut1",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # nnmut3 = circleButton(
  #   inputId = "nnmut3",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # 
  # nnexp1 = circleButton(
  #   inputId = "nnexp1",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # nnexp3 = circleButton(
  #   inputId = "nnexp3",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # 
  
  
  ####################################################################################################
  #Tab2:Prediction:
  #####Upload both files#####
  submitmodel.p = actionBttn(
    inputId = "submitmodel.p",
    label = "Submit",
    size = "md",
    icon = icon("play-circle"),
    style = "material-flat",
    color = "royal"
  ),
  submitexample.p = actionBttn(
    inputId = "submitexample.p",
    label = "Run Example",
    size = "md",
    icon = icon("rocket"),
    style = "material-flat",
    color = "primary"
  ),
  #example ID--> example.pd.both
  
  mut_upload.p = fileInput(
    "mut_data.p",
    NULL,
    multiple = FALSE,
    accept = c("text/csv",".maf",
               "text/comma-separated-values,text/plain",
               ".csv"),
    width = "100%"
  ),
  
  gene_list.p = textAreaInput(
    "gene_list.p",
    label = "",
    width = "100%",
    height = "200px",
    placeholder = "Paste gene symbols here as separate lines, for example:\nCDK4\nCDK6\nCCND1"
  ),
  exp_upload.p = fileInput(
    "exp_data.p",
    NULL,
    multiple = FALSE,
    accept = c("text/csv",
               "text/comma-separated-values,text/plain",
               ".csv"),
    width = "100%"
  ),
  
  exp_type.p = prettyRadioButtons(
    "exp_type.p",
    label = "Data Type",
    choices =
      c("FPKM", "TPM"),
    inline = TRUE,
    shape = "square",
    selected =  "TPM"
  ),
  
  log_transform.p = prettyRadioButtons(
    "log_transform.p",
    label = "Log Transformed (log2(data+1))?",
    choiceNames = c("Yes", "No"),
    choiceValues = c(1, 2),
    inline = TRUE,
    shape = "square",
    selected = 1
  ),
  
  
  cell_type.p = prettyRadioButtons(
    "cell_type1.p",
    label = "",
    choiceNames = c("Cell line", "Tumor"),
    choiceValues = c("CCLE", "TCGA"),
    inline = TRUE,
    shape = "square",
    selected = "CCLE"
  ),
  
  
  #Show results:
  # Table2-1: pd.table.both
  # Figure2-2:pd.waterfall.both
  # Figure2-3:pd.structure.both
  # Text2-3: pd.drug.info.both
  # Figure2-4:pd.density.both
  # Figure2-5: pd.box.plot.both
  # Figure2-6: pd.network.both
  
  
  pd.table.both = shinycssloaders::withSpinner(dataTableOutput("pd.table.both")),
  pd.waterfall.both = shinycssloaders::withSpinner(plotlyOutput("pd.waterfall.both")),
  pd.structure.both = shinycssloaders::withSpinner(imageOutput("pd.structure.both")),
  pd.drug.info.both = shinycssloaders::withSpinner(plotOutput("pd.drug.info.both")),
  pd.density.both = shinycssloaders::withSpinner(plotlyOutput("pd.density.both")),
  pd.box.plot.both = shinycssloaders::withSpinner(plotlyOutput("pd.box.plot.both")),
  pd.network.both = shinycssloaders::withSpinner(visNetworkOutput("pd.network.both")),
  
  #new figure to show sub box plot
  pd.cancer.select.both = selectInput(
    inputId = "pd.cancer.select.both",
    label = "Select a Cancer Type",
    choices = character(0),
    width = "100%"
  ),
  pd.sub.box.plot.both = shinycssloaders::withSpinner(plotlyOutput("pd.sub.box.plot.both")),
  
  pd.network.both = shinycssloaders::withSpinner(visNetworkOutput("pd.network.both")),
  
  
  #####Upload mutation only#####
  submitmodel2.p = actionBttn(
    inputId = "submitmodel2.p",
    label = "Submit",
    size = "md",
    icon = icon("play-circle"),
    style = "material-flat",
    color = "royal"
  ),
  submitexample2.p = actionBttn(
    inputId = "submitexample2.p",
    label = "Run Example",
    size = "md",
    icon = icon("rocket"),
    style = "material-flat",
    color = "primary"
  ),
  #example ID -->example.pd.mut
  
  mut_upload2.p = fileInput(
    "mut_data2.p",
    NULL,
    multiple = FALSE,
    accept = c("text/csv",".maf",
               "text/comma-separated-values,text/plain",
               ".csv"),
    width = "100%"
  ),
  
  gene_list2.p = textAreaInput(
    "gene_list2.p",
    label = "",
    width = "100%",
    height = "200px",
    placeholder = "Paste gene symbols here as separate lines, for example:\nCDK4\nCDK6\nCCND1"
  ),
  cell_type2.p = prettyRadioButtons(
    "cell_type2.p",
    label = "",
    choiceNames  = c("Cell line", "Tumor"),
    choiceValues = c("CCLE", "TCGA"),
    inline = TRUE,
    shape = "square",
    selected = "CCLE"
  ),
  #Show results:
  # Table2-1: pd.table.mut
  # Figure2-2:pd.waterfall.mut
  # Figure2-3:pd.structure.mut
  # Text2-3: pd.drug.info.mut
  # Figure2-4:pd.density.mut
  # Figure2-5: pd.box.plot.mut
  # Figure2-6: pd.network.mut
  
  pd.table.mut = shinycssloaders::withSpinner(dataTableOutput("pd.table.mut")),
  pd.waterfall.mut = shinycssloaders::withSpinner(plotlyOutput("pd.waterfall.mut")),
  pd.structure.mut = shinycssloaders::withSpinner(imageOutput("pd.structure.mut")),
  pd.drug.info.mut = shinycssloaders::withSpinner(plotOutput("pd.drug.info.mut")),
  pd.density.mut = shinycssloaders::withSpinner(plotlyOutput("pd.density.mut")),
  pd.box.plot.mut = shinycssloaders::withSpinner(plotlyOutput("pd.box.plot.mut")),
  pd.network.mut = shinycssloaders::withSpinner(visNetworkOutput("pd.network.mut")),
  
  #new subplots for mutation
  pd.cancer.select.mut = selectInput(
    inputId = "pd.cancer.select.mut",
    label = "Select a Cancer Type",
    choices = character(0),
    width = "100%"
  ),
  pd.sub.box.plot.mut = shinycssloaders::withSpinner(plotlyOutput("pd.sub.box.plot.mut")),
  
  pd.network.mut = shinycssloaders::withSpinner(visNetworkOutput("pd.network.mut")),
  
  #####upload gene expression only#####
  submitmodel3.p = actionBttn(
    inputId = "submitmodel3.p",
    label = "Submit",
    size = "md",
    icon = icon("play-circle"),
    style = "material-flat",
    color = "royal"
  ),
  submitexample3.p = actionBttn(
    inputId = "submitexample3.p",
    label = "Run Example",
    size = "md",
    icon = icon("rocket"),
    style = "material-flat",
    color = "primary"
  ),
  
  #example ID -->example.pd.exp
  
  exp_upload2.p = fileInput(
    "exp_data2.p",
    NULL,
    multiple = FALSE,
    accept = c("text/csv",
               "text/comma-separated-values,text/plain",
               ".csv"),
    width = "100%"
  ),
  exp_type2.p = prettyRadioButtons(
    "exp_type2.p",
    label = "Data Type",
    choices =
      c("FPKM", "TPM"),
    inline = TRUE,
    shape = "square",
    selected = "TPM"
  ),
  log_transform2.p = prettyRadioButtons(
    "log_transform2.p",
    label = "Log Transformed (log2(data+1))?",
    choiceNames = c("Yes", "No"),
    choiceValues = c(1, 2),
    inline = TRUE,
    shape = "square",
    selected = 1
  ),
  cell_type3.p = prettyRadioButtons(
    "cell_type3.p",
    label = "",
    choiceNames  = c("Cell line", "Tumor"),
    choiceValues = c("CCLE", "TCGA"),
    inline = TRUE,
    shape = "square",
    selected = "CCLE"
  ),
  
  #Show results:
  # Table2-1: pd.table.exp
  # Figure2-2:pd.waterfall.exp
  # Figure2-3:pd.structure.exp
  # Text2-3: pd.drug.info.exp
  # Figure2-4:pd.density.exp
  # Figure2-5: pd.box.plot.exp
  # Figure2-6: pd.network.exp
  
  pd.table.exp = shinycssloaders::withSpinner(dataTableOutput("pd.table.exp")),
  pd.waterfall.exp = shinycssloaders::withSpinner(plotlyOutput("pd.waterfall.exp")),
  pd.structure.exp = shinycssloaders::withSpinner(imageOutput("pd.structure.exp")),
  pd.drug.info.exp = shinycssloaders::withSpinner(plotOutput("pd.drug.info.exp")),
  pd.density.exp = shinycssloaders::withSpinner(plotlyOutput("pd.density.exp")),
  pd.box.plot.exp = shinycssloaders::withSpinner(plotlyOutput("pd.box.plot.exp")),
  
  #sub box plot for expression
  pd.cancer.select.exp = selectInput(
    inputId = "pd.cancer.select.exp",
    label = "Select a Cancer Type",
    choices = character(0),
    width = "100%"
  ),
  pd.sub.box.plot.exp = shinycssloaders::withSpinner(plotlyOutput("pd.sub.box.plot.exp")),
  
  pd.network.exp = shinycssloaders::withSpinner(visNetworkOutput("pd.network.exp"))
  
  #helptabs (please alter these for this year's app)
  # helptab3p1 = actionLink(inputId = "helptab3p1", label = "Input Data Format"),
  # helptab1p1 = actionLink(inputId = "helptab1p1", label = "Find Similar Samples (Find Sample)"),
  # helptab2p1 = actionLink(inputId = "helptab2p1", label = "Predict Drug Response (Find Drug)"),
  # 
  # helptab3p2 = actionLink(inputId = "helptab3p2", label = "Input Data Format"),
  # helptab1p2 = actionLink(inputId = "helptab1p2", label = "Find Similar Samples (Find Sample)"),
  # helptab2p2 = actionLink(inputId = "helptab2p2", label = "Predict Drug Response (Find Drug)"),
  # 
  # helptab3p3 = actionLink(inputId = "helptab3p3", label = "Input Data Format"),
  # helptab1p3 = actionLink(inputId = "helptab1p3", label = "Find Similar Samples (Find Sample)"),
  # helptab2p3 = actionLink(inputId = "helptab2p3", label = "Predict Drug Response (Find Drug)"),
  # 
  # helptab3p4 = actionLink(inputId = "helptab3p4", label = "Input Data Format"),
  # helptab1p4 = actionLink(inputId = "helptab1p4", label = "Find Similar Samples (Find Sample)"),
  # helptab2p4 = actionLink(inputId = "helptab2p4", label = "Predict Drug Response (Find Drug)"),
  # 
  # helptab3p5 = actionLink(inputId = "helptab3p5", label = "Input Data Format"),
  # helptab1p5 = actionLink(inputId = "helptab1p5", label = "Find Similar Samples (Find Sample)"),
  # helptab2p5 = actionLink(inputId = "helptab2p5", label = "Predict Drug Response (Find Drug)"),
  # 
  # help_input1 = actionLink(
  #   inputId = "help_input1",
  #   label = "Learn more about input data format and download examples",
  #   style = "color:blueviolet;
  #                                        border-bottom: 1px solid #a9529e;"
  # ),
  # help_input2 = actionLink(
  #   inputId = "help_input2",
  #   label = "Learn more about input data format and download examples",
  #   style = "color:blueviolet;
  #                                        border-bottom: 1px solid #a9529e;"
  # ),
  # 
  # pdboth1 = circleButton(
  #   inputId = "pdboth1",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # pdboth2 = circleButton(
  #   inputId = "pdboth2",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # pdboth3 = circleButton(
  #   inputId = "pdboth3",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # pdboth5 = circleButton(
  #   inputId = "pdboth5",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # pdboth6 = circleButton(
  #   inputId = "pdboth6",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # 
  # #NEW MODAL
  # pdboth6sub = circleButton(
  #   inputId = "pdboth6sub",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # pdboth7 = circleButton(
  #   inputId = "pdboth7",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # 
  # pdmut1 = circleButton(
  #   inputId = "pdmut1",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # pdmut2 = circleButton(
  #   inputId = "pdmut2",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # pdmut3 = circleButton(
  #   inputId = "pdmut3",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # pdmut5 = circleButton(
  #   inputId = "pdmut5",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # pdmut6 = circleButton(
  #   inputId = "pdmut6",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # #NEW MODAL
  # pdmut6sub = circleButton(
  #   inputId = "pdmut6sub",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # pdmut7 = circleButton(
  #   inputId = "pdmut7",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # 
  # pdexp1 = circleButton(
  #   inputId = "pdexp1",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # pdexp2 = circleButton(
  #   inputId = "pdexp2",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # pdexp3 = circleButton(
  #   inputId = "pdexp3",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # pdexp5 = circleButton(
  #   inputId = "pdexp5",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # pdexp6 = circleButton(
  #   inputId = "pdexp6",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # #NEW MODAL
  # pdexp6sub = circleButton(
  #   inputId = "pdexp6sub",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # ),
  # pdexp7 = circleButton(
  #   inputId = "pdexp7",
  #   icon = icon("question"),
  #   style = "color:darkblue; margin-bottom: 0.9vh;"
  # )
  
)
