# setwd("~/Desktop/DeepDRShiny_version2/")
# 
# inputTable.ori <- readRDS("data/NN_table_tcga.RDS")
# inputTable <- inputTable.ori[1:10,]
# target = "User_Sample"
# 
# cell_type = "CCLE"
# 
# gdsc_ic50 = readRDS("data/GDSC_paper/gdsc_and_ccle_704_overlapped_cells_ic50.RDS")
# 
# tcga.predict50 <- readRDS("data/TCGA/tcga_predict_ic50_9059samples.RDS")

ShinyDeepDR_visNetwork_nn <- function(inputTable, 
                                      target,
                                      cell_type = c("CCLE", "TCGA"),
                                      gdsc_ic50, tcga.predict50, 
                                      num_drug = 10){
  library("visNetwork")
  library(dplyr)
  library(tibble)
  
  if(cell_type == "CCLE"){
    rownames(gdsc_ic50) <- gdsc_ic50$DepMapID
    gdsc_ic50.t <- t(gdsc_ic50[,-c(1,2)])
    
    #First layer
    shape_1 <- c("dot",rep("dot",nrow(inputTable)))
    color_1 <- c("blue", rep( "peachpuff", nrow(inputTable)))
    
    nodes <- data.frame(id = seq.int(from =1,to = (nrow(inputTable)+1), by = 1),
                        label = c(target, inputTable[,"cell_line_name"]), 
                        cell_id = c(target, inputTable[,"DepMapID"]),
                        shape = shape_1,
                        color = color_1, stringsAsFactors = F)
    
    if(sum(colnames(inputTable) %in% "Corr") == 1){
      edges <- data.frame( id = seq.int(from = 1,to = nrow(inputTable), by = 1),
                           from = rep(1,nrow(inputTable)), 
                           to = seq.int(from = 2, to = (nrow(inputTable)+1), by = 1),
                           width = 10*(inputTable[, 2]),
                           length = rep(200,nrow(inputTable)),
                           color = rep('black',nrow(inputTable)),
                           stringsAsFactors = F)
    }else{
      edges <- data.frame( id = seq.int(from = 1,to = nrow(inputTable), by = 1),
                           from = rep(1,nrow(inputTable)), 
                           to = seq.int(from = 2, to = (nrow(inputTable)+1), by = 1),
                           width = 50*(inputTable[, 2]),
                           length = rep(200,nrow(inputTable)),
                           color = rep('black',nrow(inputTable)),
                           stringsAsFactors = F)
    }


    
    if(sum(edges$width < 0) >0){
      edges$color[which(edges$width < 0)] <- "dodgerblue"
    }
  
    #2nd layer: select top 10 drug with small IC50 for each samples
    for (i in 2:nrow(nodes)) {
      tb <- data.frame(t(gdsc_ic50[which(gdsc_ic50$DepMapID == nodes$cell_id[i]),-c(1,2)]))
      tb <- tb %>% add_column(DrugName = rownames(tb),.before = colnames(tb)[1])
      #tb <- tb[which(!tb[,2] == "NaN"),]
      tb.sort <- as.data.frame(tb[order(tb[,2]),])
      
      if(i==2){
        top10Drug <- tb.sort$DrugName[1:num_drug]
       # cat("Sample: ", nodes$cell_id[i],"\ drug: ", paste(top10Drug ,collapse = ","), "\n" )
      }else{
        top10Drug.sub <-tb.sort$DrugName[1:num_drug]
        top10Drug <- c(top10Drug,top10Drug.sub)
       # cat("Sample: ", nodes$cell_id[i],"\ drug: ", paste(top10Drug.sub ,collapse = ","), "\n" )
      }
      
      top10Drug.out <- unique(top10Drug)

    }
    
    nodes.final <- rbind(nodes,
                   data.frame(id = seq.int(from=(nrow(nodes)+1),to =( nrow(nodes)+length(top10Drug.out)), by = 1),
                              label = top10Drug.out,
                              cell_id = NA,
                              shape = rep("diamond", length(top10Drug.out)),
                              color = rep("skyblue", length(top10Drug.out)),stringsAsFactors = F))
    
    for (i in 2:nrow(nodes)) {
      tb <- data.frame(t(gdsc_ic50[which(gdsc_ic50$DepMapID == nodes$cell_id[i]),-c(1,2)]))
      tb <- tb %>% add_column(DrugName = rownames(tb),.before = colnames(tb)[1])
      #tb <- tb[which(!tb[,2] == "NaN"),]
      tb.sort <- as.data.frame(tb[order(tb[,2]),])
      
      if(i==2){
        top10Drug <- tb.sort$DrugName[1:num_drug]
        edge.sub <- data.frame(matrix(data = NA, nrow = length(top10Drug), ncol = ncol(edges)), stringsAsFactors = F)
        colnames(edge.sub) <- colnames(edges)
        edge.sub$id <- seq.int(from=(nrow(edges)+1), to = (nrow(edges)+length(top10Drug)),by = 1)
        edge.sub$from <- i
        edge.sub$length <- 200
        edge.sub$color <- "lightcoral"
          
        for (j in 1:length(top10Drug)) {
         edge.sub$to[j] <- nodes.final$id[which(nodes.final$label == top10Drug[j])]
         edge.sub$width[j] <- -tb.sort[top10Drug[j],2]
          
          
        }
        
        edges.final <- rbind(edges, edge.sub)
        #cat("Sample: ", nodes$cell_id[i],"\ drug: ", paste(top10Drug ,collapse = ","), "\n" )
      }else{
        
        top10Drug <- tb.sort$DrugName[1:num_drug]
        edge.sub <- data.frame(matrix(data = NA, nrow = length(top10Drug), ncol = ncol(edges.final)), stringsAsFactors = F)
        colnames(edge.sub) <- colnames(edges.final)
        edge.sub$id <- seq.int(from=(nrow(edges)+1), to = (nrow(edges)+length(top10Drug)),by = 1)
        edge.sub$from <- i
        edge.sub$length <- 200
        edge.sub$color <- "lightcoral"
        
        for (j in 1:length(top10Drug)) {
          edge.sub$to[j] <- nodes.final$id[which(nodes.final$label == top10Drug[j])]
          edge.sub$width[j] <- -tb.sort[top10Drug[j],2]
          
          
        }
       # cat("Sample: ", nodes$cell_id[i],"\ drug: ", paste(top10Drug.sub ,collapse = ","), "\n" )
      }
      edges.final <- rbind(edges.final, edge.sub)
      
    }
    
  }
  
  if(cell_type == "TCGA"){
    


    #First layer
    shape_1 <- c("dot",rep("dot",nrow(inputTable)))
    color_1 <- c("blue", rep( "peachpuff", nrow(inputTable)))
    
    nodes <- data.frame(id = seq.int(from =1,to = (nrow(inputTable)+1), by = 1),
                        label = c(target, inputTable[,"PatientID" ]), 
                        cell_id = c(target, inputTable[,"PatientID" ]),
                        shape = shape_1,
                        color = color_1)
    
    if(sum(colnames(inputTable) %in% "Corr") == 1){
      edges <- data.frame( id = seq.int(from = 1, to = nrow(inputTable), by = 1),
                           from = rep(1,nrow(inputTable)), 
                           to = seq.int(from = 2, to = (nrow(inputTable)+1), by = 1),
                           width = 5*(inputTable[, 2]),
                           length = rep(200,nrow(inputTable)),
                           color = rep('black',nrow(inputTable)))
    }else{
      edges <- data.frame( id = seq.int(from = 1, to = nrow(inputTable), by = 1),
                           from = rep(1,nrow(inputTable)), 
                           to = seq.int(from = 2, to = (nrow(inputTable)+1), by = 1),
                           width = 20*(inputTable[, 2]),
                           length = rep(200,nrow(inputTable)),
                           color = rep('black',nrow(inputTable)))
    }

    
    if(sum(edges$width < 0) >0){
      edges$color[which(edges$width < 0)] <- "dodgerblue"
    }
    
    #2nd layer: select top 10 drug with small IC50 for each samples
    for (i in 2:nrow(nodes)) {
      tb <- data.frame(tcga.predict50[,c(1,which(substr(colnames(tcga.predict50),1,12 )== nodes$cell_id[i]))])
      tb.sort <- as.data.frame(tb[order(tb[,2]),])
      
      if(i==2){
        top10Drug <- tb.sort$drug_names[1:num_drug]
        #cat("Sample: ", nodes$cell_id[i],"\ drug: ", paste(top10Drug ,collapse = ","), "\n" )
      }else{
        top10Drug.sub <-tb.sort$drug_names[1:num_drug]
        top10Drug <- c(top10Drug,top10Drug.sub)
        #cat("Sample: ", nodes$cell_id[i],"\ drug: ", paste(top10Drug.sub ,collapse = ","), "\n" )
      }
      
      top10Drug.out <- unique(top10Drug)
      
    }
    
    nodes.final <- rbind(nodes,
                         data.frame(id = seq.int(from = (nrow(nodes)+1),to = (nrow(nodes)+length(top10Drug.out)), by = 1),
                                    label = top10Drug.out,
                                    cell_id = NA,
                                    shape = rep("diamond", length(top10Drug.out)),
                                    color = rep("skyblue", length(top10Drug.out)),stringsAsFactors = F))
    
    for (i in 2:nrow(nodes)) {
      tb <- data.frame(tcga.predict50[,c(1,which(substr(colnames(tcga.predict50),1,12 )== nodes$cell_id[i]))])
      tb.sort <- as.data.frame(tb[order(tb[,2]),])
      
      if(i==2){
        top10Drug <-tb.sort$drug_names[1:num_drug]
        edge.sub <- data.frame(matrix(data = NA, nrow = length(top10Drug), ncol = ncol(edges)), stringsAsFactors = F)
        colnames(edge.sub) <- colnames(edges)
        edge.sub$id <- seq.int(from = (nrow(edges)+1),to = ( nrow(edges)+length(top10Drug)),by = 1)
        edge.sub$from <- i
        edge.sub$length <- 200
        edge.sub$color <- "lightcoral"
        
        for (j in 1:length(top10Drug)) {
          edge.sub$to[j] <- nodes.final$id[which(nodes.final$label == top10Drug[j])]
          edge.sub$width[j] <- -tb.sort[which(tb.sort$drug_names == top10Drug[j]),2]
          
          
        }
        
        edges.final <- rbind(edges, edge.sub)
        cat("Sample: ", nodes$cell_id[i],"\ drug: ", paste(top10Drug ,collapse = ","), "\n" )
      }else{
        
        top10Drug <- tb.sort$drug_names[1:num_drug]
        edge.sub <- data.frame(matrix(data = NA, nrow = length(top10Drug), ncol = ncol(edges.final)), stringsAsFactors = F)
        colnames(edge.sub) <- colnames(edges.final)
        edge.sub$id <- seq.int(from = (nrow(edges.final)+1),to = ( nrow(edges.final)+length(top10Drug)),by = 1)
        edge.sub$from <- i
        edge.sub$length <- 200
        edge.sub$color <- "lightcoral"
        
        for (j in 1:length(top10Drug)) {
          edge.sub$to[j] <- nodes.final$id[which(nodes.final$label == top10Drug[j])]
          edge.sub$width[j] <-  -tb.sort[which(tb.sort$drug_names == top10Drug[j]),2]
          
          
        }
        #cat("Sample: ", nodes$cell_id[i],"\ drug: ", paste(top10Drug.sub ,collapse = ","), "\n" )
      }
      edges.final <- rbind(edges.final, edge.sub)
      
    }
    
    
    
  }
  
  
  if(sum(colnames(inputTable) %in% "Corr") == 1){
    label = c("Top similarity (correlation)","Top sensitivity")
  }else{
    label = c("Top similarity (Jaccard index)","Top sensitivity")
  }
  
  edges.final$id <- seq.int(from = 1,to =  nrow(edges.final), by = 1)
  nodes.final$group <- c("Query sample",rep("Similar sample", nrow(nodes)-1) ,rep("Drug",length(top10Drug.out)))
  ledges <- data.frame(color = c("black", "lightcoral"), 
                       label = label,
                       arrows= list(to = list(enabled = F, 
                                              scaleFactor = 2, type = 'l'))
                       ) #Remove arrow of legend
 
  p <- visNetwork::visNetwork(nodes = nodes.final,edges = edges.final)%>%
    visIgraphLayout(layout = "layout_with_graphopt") %>%
    visLayout(randomSeed = 12)%>%
    visEdges(smooth = FALSE) %>%
    visNodes(size=40, shadow = list(enabled = TRUE, size = 10)) %>%
    visGroups(groupname = "Query sample" , color = "blue")%>%
    visGroups(groupname = "Similar sample" , color = "peachpuff")%>%
    visGroups(groupname = "Drug" , color = "skyblue")%>%
      visOptions(selectedBy = "group", 
                 highlightNearest = TRUE, 
                 nodesIdSelection = TRUE,
                 autoResize = T) %>%
      visPhysics(stabilization = FALSE)%>%
    visLegend(width=0.25, position = "right", useGroups = F, addEdges = ledges, 
              addNodes = list(
                list(label = "Query sample", shape = "dot", color = "blue"),
                list(label = "Similar sample", shape = "dot", color = "peachpuff"),
                list(label = "Drug", shape = "diamond", color = "skyblue")
                ))


  return(p)
  
}
