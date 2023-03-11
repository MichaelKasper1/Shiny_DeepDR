# Table2-1: Table
ShinyDeepDR_DrugResponsePercentile<- function (drug_scores){
  last_element = length(drug_scores)
  reference_drugScores = drug_scores[1:(last_element-1)]
  reference_drugScores = reference_drugScores[!is.na(reference_drugScores)]
  user_drugScore = drug_scores[last_element]
  p = length(reference_drugScores[reference_drugScores <= user_drugScore])/length(reference_drugScores)*100
  return(p)
}

ShinyDeepDR_predict_drug_table <- function(inputData,
                                           cell_type,
                                           ccle_ic50Pred,
                                           gdsc_ic50,
                                           tcga_ic50Pred,
                                           gdsc_drugInfo) {
  if (cell_type == "TCGA"){
    
    real_ic50_range = t(apply(tcga_ic50Pred[, 2:9060], 1, range, na.rm=TRUE))
    #real_ic50_range = formatC(real_ic50_range, format = "e", digits = 2)
    
    temp = cbind(tcga_ic50Pred[,2:9060], inputData$Log_predict_IC50)
    real_ic50_percentile = apply(temp, 1, ShinyDeepDR_DrugResponsePercentile)
    
    drugPredTable = cbind(inputData, real_ic50_range, round(real_ic50_percentile,2))
    # drugPredTable$Log_predict_IC50 = formatC(drugPredTable$Log_predict_IC50, format = "e", digits = 2)
    colnames(drugPredTable) = c('DrugName', 'Log_predict_IC50', 'TCGA-PredIC50 min',
                                'TCGA-PredIC50 max', 'TCGA-PredIC50 Percentile')
    
  } else {
    gdsc_ic50trans = data.frame(t(gdsc_ic50[,3:267]))
    real_ic50_range = t(apply(gdsc_ic50trans, 1, range, na.rm=TRUE))
    #real_ic50_range = formatC(real_ic50_range, format = "e", digits = 2)
    
    temp = cbind(gdsc_ic50trans, inputData$Log_predict_IC50)
    real_ic50_percentile = apply(temp, 1, ShinyDeepDR_DrugResponsePercentile)
    
    ccle_ic50trans = data.frame(t(ccle_ic50Pred))
    ccleP_ic50_range = t(apply(ccle_ic50trans, 1, range, na.rm=TRUE))
    #ccleP_ic50_range = formatC(ccleP_ic50_range, format = "e", digits = 2)
    
    temp = cbind(ccle_ic50trans, inputData$Log_predict_IC50)
    ccleP_ic50_percentile = apply(temp, 1, ShinyDeepDR_DrugResponsePercentile)
    
    drugPredTable = cbind(inputData,
                          round(real_ic50_range, digits = 2), 
                          round(real_ic50_percentile,2),
                          round(ccleP_ic50_range, digits = 2), 
                          round(ccleP_ic50_percentile,2))
    #drugPredTable$Log_predict_IC50 = formatC(drugPredTable$Log_predict_IC50, format = "e", digits = 2)
    colnames(drugPredTable) = c('DrugName', 'Log_predict_IC50', 'GDSC-CCLE IC50 min',
                                'GDSC-CCLE IC50 max', 'GDSC-CCLE IC50 Percentile',
                                'CCLE-Pred IC50 min', 'CCLE-Pred IC50 max',
                                'CCLE-Pred IC50 Percentile')
    
  }
  drug.info.drugname <- gdsc_drugInfo$colname
  drug.info.drugname[256] <- "Salubrinal"
  drugPredTable <- drugPredTable[match(drugPredTable$DrugName, drug.info.drugname),]
  drugPredTable = cbind(drugPredTable, gdsc_drugInfo[,3:ncol(gdsc_drugInfo)])
  drugPredTable$Log_predict_IC50 <- round(drugPredTable$Log_predict_IC50, digits = 2)
  
  return(drugPredTable)
}




# Figure2-2: Waterfall plot  
ShinyDeepDR_WaterFall <- function(inputData,num_drugs = 265 ) {

  
  drugnames = make.names(inputData$DrugName, unique = TRUE)
  ll = order(inputData$Log_predict_IC50)
  ord_DN = as.matrix(drugnames[ll])
  ord_DR = as.matrix(as.numeric(inputData$Log_predict_IC50[ll]))
  y = ord_DR
  measure = c(rep('absolute', num_drugs))
  data = data.frame(x=factor(ord_DN,levels = ord_DN), y, measure)
  
  
  fig <- plot_ly(data, x = ~x, y = ~y, measure = ~measure, type = "waterfall", 
                 decreasing = list(marker = list(color = "Maroon", 
                                                 line = list(color = "red", width = 2))),
                 increasing = (marker = list(color = "Green")),
                 totals = list(marker = list(color = "blue", 
                                             line = list(color = 'skyblue', width = 3))))
  
  fig <- fig %>% layout(title = " ",
                        xaxis = list(title = ""),
                        yaxis = list(title = "Predicted log(IC50)"))
  
  return(fig)
}

# Figure 2-4:
ShinyDeepDR_vline <- function(x = 0, y_height, color="red") {
  list(type = "line", 
       y0 = 0, 
       y1 = y_height, 
       yref = "y",
       x0 = x, 
       x1 = x, 
       showlegend = TRUE,
       line = list(color = color)
  )
}

ShinyDeepDR_DensityPlot <- function(inputData, cell_type, select.drug,ccle_ic50Pred, gdsc_ic50,
                                    tcga_ic50Pred ) {
  
  match_id = which(colnames(gdsc_ic50)[-c(1,2)] == select.drug)
  
  if (cell_type == "CCLE"){
    ccle_ic50data = data.frame(t(ccle_ic50Pred))
    dataC = ccle_ic50data[match_id[1],]
    dataC = dataC[!is.na(dataC)]
    density_ccle = density(as.numeric(dataC))
    fig <- plot_ly(alpha = 0.3)
    fig <- fig %>% add_trace(x = round(density_ccle$x, 3),
                             y = round(density_ccle$y, 3),
                             type = "scatter",
                             fill = "tozeroy",
                             mode = "lines",
                             line = list(color = "#99ccff"),
                             fillcolor = 'rgba(153, 206, 235, 0.5)',
                             name = "DeepDR prediction (622 cell lines)")
    
    gdsc_ic50data = data.frame(t(gdsc_ic50[,3:267]))
    
    data = gdsc_ic50data[match_id[1],]
    data = data[!is.na(data)]
    pre_ic50 = as.numeric(inputData$Log_predict_IC50[which(inputData$DrugName == select.drug)])
    density_gdsc = density(as.numeric(data))
    fig <- fig %>% add_trace(x = round(density_gdsc$x, 3),
                             y = round(density_gdsc$y, 3),
                             type = "scatter",
                             fill = "tozeroy",
                             mode = "lines",
                             line = list(color = "#FFDAB9"),
                             fillcolor = 'rgba(255, 218, 185, 0.5)',
                             name = "GDSC data (622 cell lines)")
    
    
    
    
    y_max = max(density_gdsc$y, density_ccle$y)
    fig <- fig %>% layout(showlegend=TRUE, legend=list(title=list(text='Reference')))
    fig <- fig %>% layout(shapes = list(ShinyDeepDR_vline(x = pre_ic50, y_height = 0.5*y_max)))
    fig <- fig %>% layout(xaxis=list(title="Predicted log(IC50)", titlefont=list(size=15)),
                          yaxis=list(title="Density of samples", titlefont=list(size=15)))
    fig <- fig %>% add_text(showlegend=FALSE, x=pre_ic50, y=0.6*y_max, 
                            name=select.drug,
                            text=select.drug,  hovertemplate=paste(select.drug))
    
  } else {
    
    data = tcga_ic50Pred[match_id[1], 2:9060]
    data = data[!is.na(data)]
    pre_ic50 = as.numeric(inputData$Log_predict_IC50[which(inputData$DrugName == select.drug)])
    density_tcga = density(as.numeric(data))
    
    fig <- plot_ly(alpha = 0.3)
    fig <- fig %>% add_trace(x = round(density_tcga$x, 3),
                             y = round(density_tcga$y, 3),
                             type = "scatter",
                             fill = "tozeroy",
                             mode = "lines",
                             line = list(color = "#BEBADA"),
                             fillcolor = 'rgba(190, 186, 218, 0.5)',
                             name = "DeepDR prediction (9059 tumors)")
    
    y_max = max(density_tcga$y)
    fig <- fig %>% layout(showlegend=TRUE, legend=list(title=list(text='Reference')))
    fig <- fig %>% layout(shapes = list(ShinyDeepDR_vline(x = pre_ic50,y_height =  0.5*y_max)))
    fig <- fig %>% layout(xaxis=list(title="Predicted log(IC50)", titlefont=list(size=15)),
                          yaxis=list(title="Density of samples", titlefont=list(size=15)))
    fig <- fig %>% add_text(showlegend=FALSE, x=pre_ic50, y=0.6*y_max, 
                            name=select.drug,
                            text=select.drug,  hovertemplate=paste(select.drug))
  }
  return(fig)
}


# Figure 2-5: PanCan BoxPlot
ShinyDeepDR_hline <- function(y , x_width, color="blue") {
  list(type = "line", 
       y0 = y, 
       y1 = y, 
       xref = "x",
       x0 = 0, 
       x1 = x_width, 
       line = list(color = color, dash="dot")
  )
}

# inputData = table.out
# select.drug = "Docetaxel"
# ccle_ic50Pred = ccle_ic50_Pred[,table.out$DrugName]
# gdsc_ic50 = drug.ic50
# gdsc_cellInfo = gdsc_ccle_cellInfo
  
ShinyDeepDR_BoxPlot <- function(inputData, cell_type, select.drug, ccle_ic50Pred, gdsc_ic50, gdsc_cellInfo,
                                tcga_ic50Pred, tcga_clinicalData){
  
  match_id = match(select.drug, inputData$DrugName)
  pre_ic50 = as.numeric(inputData$Log_predict_IC50[match_id[1]])
  
  if (cell_type == 'CCLE'){
    colnames(gdsc_ic50)[258] <- "Salubrinal"
    # data_temp = as.matrix(t(gdsc_ic50[match_id[1],]))
    data_temp = as.matrix(t(gdsc_ic50[,select.drug]))
    data_gdsc = data.frame('PanCan'= gdsc_cellInfo$primary_disease, 'logIC50'= data_temp[1,], 
                           'source' = replicate(length(gdsc_cellInfo$primary_disease),'GDSC data'))
    
    data_temp = as.matrix(ccle_ic50Pred[,match_id[1]])
    data_ccle = data.frame('PanCan'= gdsc_cellInfo$primary_disease, 'logIC50'= data_temp[,1], 
                           'source' = replicate(length(gdsc_cellInfo$primary_disease),'DeepDR prediction'))
    
    unique_tumors = unique(gdsc_cellInfo$primary_disease)
    num_unique_tumor = length(unique_tumors)
    
    data_merge = rbind(data_gdsc, data_ccle)
    
    cc = as.matrix(table(gdsc_cellInfo$primary_disease))
    x_text = apply(t(rownames(cc)), 1, paste,"(",cc[,1],")")
    
    fig = plot_ly(alpha = 0.2, data_merge, x=~PanCan, y=~logIC50, color=~source, 
                  type="box", colors = c("skyblue", 'peachpuff'))
    
    fig = fig %>% layout(boxmode="group", xaxis=list(title="", 
                                                     categoryarray='array',
                                                     categoryorder="category ascending",
                                                     showticklabels=TRUE, 
                                                     tickmode='array', 
                                                     ticktext=x_text,
                                                     tickvals=sort(unique_tumors)),
                         yaxis=list(title = "Predicted log(IC50)", titlefont=list(size=15)))
    
    fig <- fig %>% layout(xaxis=list(showticklabels = TRUE, tickmode = 'array', 
                                     ticktext=x_text,
                                     tickvals=sort(unique_tumors)))
    
    fig <- fig %>% layout(shapes = list(ShinyDeepDR_hline(y = pre_ic50,
                                                          x_width = num_unique_tumor+1)))
    fig <- fig %>% add_text(showlegend = FALSE, x = 1+0.5, y = pre_ic50+3, 
                            name = select.drug, textfont = list(color = 'blue'),
                            text=select.drug, 
                            hovertemplate= paste(select.drug))
    
    return(fig)
    
  } else{
    
    data_temp = as.matrix(t(tcga_ic50Pred[match_id[1],2:9060]))
    data_tcga = data.frame('PanCan'= tcga_clinicalData$type, 'logIC50'= data_temp[,1], 
                           'source' = replicate(length(tcga_clinicalData$type),'DeepDR prediction'))
    
    unique_tumors = unique(tcga_clinicalData$type)
    num_unique_tumor = length(unique_tumors)
    
    cc = as.matrix(table(tcga_clinicalData$type))
    x_text = apply(t(rownames(cc)), 1, paste,"(",cc[,1],")")
    
    fig = plot_ly(alpha = 0.2, data_tcga, x=~PanCan, y=~logIC50, color=~source,
                  type="box", colors = c("pink"))
    
    fig = fig %>% layout(xaxis=list(title="", 
                                    categoryarray='array',
                                    categoryorder="category ascending",
                                    showticklabels=TRUE, 
                                    tickmode='array', 
                                    ticktext=x_text,
                                    tickvals=sort(unique_tumors)),
                         yaxis=list(title = "Predicted log(IC50)", titlefont=list(size=15)))
    
    fig <- fig %>% layout(shapes = list(ShinyDeepDR_hline(pre_ic50, num_unique_tumor)))
    fig <- fig %>% add_text(showlegend = FALSE, x = 1, y = pre_ic50+1, 
                            name = select.drug, textfont = list(color = 'blue'),
                            text=select.drug,  hovertemplate= paste(select.drug))
    return(fig)
  }
} 

# Figure 2-X: CCLE Subtype PanCan BoxPlot
shinyDeepDR_CCLESubTypeBoxPlot <- function(inputData, select.drug,
                                           pancan_type, ccle_ic50Pred,
                                           gdsc_ic50, gdsc_cellInfo){
  colnames(gdsc_ic50)[258] <- "Salubrinal"
  match_id = match(select.drug, inputData$DrugName)
  pre_ic50 = as.numeric(inputData$Log_predict_IC50[match_id[1]])
  pancan_idx = which(gdsc_cellInfo$primary_disease == pancan_type)
  
  data_temp = as.matrix(t(gdsc_ic50[match_id[1],]))
  data_temp = as.matrix(t(gdsc_ic50[,select.drug]))
  data_gdsc = data.frame('PanCan'= gdsc_cellInfo$primary_disease[pancan_idx], 
                         'Subtype' = gdsc_cellInfo$Subtype[pancan_idx],
                         'logIC50'= data_temp[1,pancan_idx], 
                         'source' = replicate(length(pancan_idx),'GDSC data'))
  
  data_temp = as.matrix(ccle_ic50Pred[,match_id[1]])
  data_ccle = data.frame('PanCan'= gdsc_cellInfo$primary_disease[pancan_idx], 
                         'Subtype' = gdsc_cellInfo$Subtype[pancan_idx],
                         'logIC50'= data_temp[pancan_idx,1], 
                         'source' = replicate(length(pancan_idx),'DeepDR prediction'))
  
  data_merge = rbind(data_gdsc, data_ccle)
  
  idx1 = which(data_merge$Subtype %in% c("", NA))
  if (!is_empty(idx1)){
    data_merge$Subtype[idx1] = "Unknown"
  }
  
  remove_idx = which(is.na(data_merge$logIC50))
  if (!is_empty(remove_idx)){
    data_merge = data_merge[-remove_idx,]
  }
  
  subtype <- data_merge$Subtype
  cc = as.matrix(table(subtype))
  x_text = apply(t(rownames(cc)), 1, paste,"(",cc[,1],")")
  
  fig = plot_ly(alpha = 0.2, data_merge, x=~Subtype, y=~logIC50, color=~source, 
                type="box", colors = c("skyblue", 'peachpuff'))
  
  fig = fig %>% layout(boxmode="group", xaxis=list(title="", 
                                                   categoryarray='array',
                                                   categoryorder="category ascending",
                                                   showticklabels=TRUE, 
                                                   tickmode='array', 
                                                   ticktext=x_text,
                                                   tickvals=sort(unique(data_merge$Subtype))),
                       yaxis=list(title = "Predicted log(IC50)",
                                  titlefont=list(size=15)))
  fig <- fig %>% layout(xaxis=list(showticklabels = TRUE, tickmode = 'array', 
                                   ticktext=x_text,
                                   tickvals=sort(unique(data_merge$Subtype))))
  fig <- fig %>% layout(shapes = list(ShinyDeepDR_hline(pre_ic50, length(unique(data_merge$Subtype))+1)))
  fig <- fig %>% add_text(showlegend = FALSE, x=1.5, y = pre_ic50+1, 
                          name = select.drug, textfont = list(color = 'blue'),
                          text=select.drug)
  return(fig)
}

# Figure 2-6: Network
ShinyDeepDR_network_for_prediction <- function(select.drug, cell_type,gdsc_ic50,
                                               tcga_ic50Pred, tcga_clinicalData,
                                               num_nodes = 10){
  
  if (cell_type == "CCLE") {
    
    drug_ic50 = gdsc_ic50[,select.drug]
    order_idx = order(drug_ic50)
    label_1 = c(select.drug, gdsc_ic50[order_idx[1:num_nodes], 1])
    color_1 = c("blue", rep( "peachpuff", num_nodes))
    shape_1 = c("diamond", rep("dot",num_nodes))
    id_1 = 1:(num_nodes+1)
    nodes <- data.frame( id = id_1, label = label_1, shape = shape_1, color = color_1 )
    
    to_1 = c(2:(num_nodes+1))
    from_1 = c(rep(1, num_nodes))
    ic50_val = gdsc_ic50[order_idx[1:num_nodes], select.drug]
    length = c( rep(150,num_nodes) )
    edges <- data.frame( id = 1:num_nodes, from = from_1, to = to_1, width = 0.5*abs(ic50_val), 
                         length=length )
    edges <- mutate(edges, color =  ifelse( ic50_val > 0,  "lightcoral","dodgerblue"))
    
    # Sub network from query drug
    cell_idx = match(nodes$label[2:(num_nodes+1)], gdsc_ic50[,1])
    
    for (n in 1:num_nodes){
      
      drug_ic50 = gdsc_ic50[cell_idx[n], 3:267]
      drug_idx = order(as.numeric(drug_ic50))
      drug_min = colnames(drug_ic50)[drug_idx[2:num_nodes]]
      drug_min_in = setdiff(drug_min,nodes$label)
      
      label_1 = c(label_1, drug_min_in)
      shape_1 = c(shape_1, rep("diamond", length(drug_min_in)))
      color_1 = c(color_1, rep("skyblue", length(drug_min_in)))
      nodes <- data.frame( id = 1:length(label_1), label = label_1, shape = shape_1,
                           color = color_1 )
      
      from_1 = c(from_1, rep(n+1, length(drug_min)))
      to_1 = c(to_1, match(drug_min, nodes$label))
      ic50_val = c(ic50_val, as.numeric(drug_ic50[drug_idx[2:num_nodes]]))
      
      # length = c( rep(200, length(to_1)) )
      length = c(length, rep(300, length(drug_min)))
      edges <- data.frame( id = 1:length(from_1), from = from_1, to = to_1, width = 0.5*abs(ic50_val), 
                           length=length )
      edges <- mutate(edges, color =  ifelse( ic50_val > 0, "lightcoral","dodgerblue"))
      
    }
    
    ledges <- data.frame(color = c("lightcoral","dodgerblue"), 
                         label = c("Top GDSC resistance (positive IC50)","Top GDSC sensitivity (negative IC50)"), 
                         arrows= list(to = list(enabled = F, 
                                                scaleFactor = 2, type = 'l')))
    nodes$group <- NA
    nodes$group[which(nodes$color == "blue")] <- "Selected drug"
    nodes$group[which(nodes$color == "skyblue")] <- "Drug"
    nodes$group[which(nodes$color == "peachpuff")] <- "Cell line"
    
    # Display the network  
    network <- visNetwork(nodes, edges) %>%
      visIgraphLayout(layout = "layout_with_graphopt") %>%
      visLayout(randomSeed = 12)%>%
      visEdges(smooth = FALSE) %>%
      visNodes(size=40, shadow = list(enabled = TRUE, size = 10) ) %>%
      visConfigure(enabled = FALSE) %>%
      visGroups(groupname = "Drug" , color = "skyblue")%>%
      visGroups(groupname = "Selected drug" , color = "blue")%>%
      visGroups(groupname = "Cell line" , color = "peachpuff")%>%
      visOptions(selectedBy = "group", 
                 highlightNearest = TRUE, 
                 nodesIdSelection = TRUE,
                 autoResize = T) %>%
      
      visLegend(width=0.25, position = "right", useGroups = F, addEdges = ledges, 
                addNodes = list(
                  list(label = "Selected drug", shape = "diamond", color = "blue"),
                  list(label = "Drug", shape = "diamond", color = "skyblue"),
                  list(label = "Cell line", shape = "dot", color = "peachpuff")))
    
    return(network)
    
  } else {
    
    drug_index = match(select.drug, tcga_ic50Pred[,1])
    # First network from query drug
    drug_ic50 = tcga_ic50Pred[drug_index, 2:9060]
    order_idx = order(as.numeric(drug_ic50))
    label_1 = c(select.drug, tcga_clinicalData$bcr_patient_barcode[order_idx[1:num_nodes]])
    color_1 = c("blue", rep( "peachpuff", num_nodes))
    shape_1 = c("diamond", rep("dot", num_nodes))
    id_1 = 1:(num_nodes+1)
    nodes <- data.frame( id = id_1, label = label_1, shape = shape_1, color = color_1 )
    
    to_1 = c(2:(num_nodes+1))
    from_1 = c(rep(1, num_nodes))
    ic50_val = as.numeric(tcga_ic50Pred[drug_index, order_idx[1:num_nodes]])
    length = c( rep(150,num_nodes) )
    edges <- data.frame( id = 1:num_nodes, from = from_1, to = to_1, width = 0.5*abs(ic50_val), 
                         length=length )
    edges <- mutate(edges, color =  ifelse( ic50_val > 0,  "lightcoral","dodgerblue"))
    
    # Sub network from query drug
    cell_idx = match(nodes$label[2:(num_nodes+1)], tcga_clinicalData$bcr_patient_barcode)
    
    for (n in 1:num_nodes){
      
      drug_ic50 = tcga_ic50Pred[, cell_idx[n]]
      drug_idx = order(as.numeric(drug_ic50))
      drug_min = tcga_ic50Pred[drug_idx[2:num_nodes],1]
      drug_min_in = setdiff(drug_min,nodes$label)
      
      label_1 = c(label_1, drug_min_in)
      shape_1 = c(shape_1, rep("diamond", length(drug_min_in)))
      color_1 = c(color_1, rep("skyblue", length(drug_min_in)))
      nodes <- data.frame( id = 1:length(label_1), label = label_1, shape = shape_1,
                           color = color_1 )
      
      from_1 = c(from_1, rep(n+1, length(drug_min)))
      to_1 = c(to_1, match(drug_min, nodes$label))
      ic50_val = c(ic50_val, as.numeric(drug_ic50[drug_idx[2:num_nodes]]))
      
      # length = c( rep(200, length(to_1)) )
      length = c(length, rep(300, length(drug_min)))
      edges <- data.frame( id = 1:length(from_1), from = from_1, to = to_1, width = 0.5*abs(ic50_val), 
                           length=length )
      edges <- mutate(edges, color =  ifelse( ic50_val > 0, "lightcoral","dodgerblue"))
      
    }
    
    ledges <- data.frame(color = c( "lightcoral", "dodgerblue"), 
                         label = c("Top predicted resistance (positive IC50)","Top predicted sensitivity (negative IC50)"), 
                         arrows= list(to = list(enabled = F, 
                                                scaleFactor = 2, type = 'l')))
    
    nodes$group <- NA
    nodes$group[which(nodes$color == "blue")] <- "Selected drug"
    nodes$group[which(nodes$color == "skyblue")] <- "Drug"
    nodes$group[which(nodes$color == "peachpuff")] <- "Tumor"
    
    # Display the network 
    network <- visNetwork(nodes, edges) %>%
      visIgraphLayout(layout = "layout_with_graphopt") %>%
      visLayout(randomSeed = 12)%>%
      visEdges(smooth = FALSE) %>%
      visNodes(size=40, shadow = list(enabled = TRUE, size = 10) ) %>%
      visConfigure(enabled = FALSE) %>%
      visGroups(groupname = "Drug" , color = "skyblue")%>%
      visGroups(groupname = "Selected drug" , color = "blue")%>%
      visGroups(groupname = "Tumor" , color = "peachpuff")%>%
      visOptions(selectedBy = "group", 
                 highlightNearest = TRUE, 
                 nodesIdSelection = TRUE,
                 autoResize = T) %>%
      visLegend(width=0.25, position = "right", useGroups = F, addEdges = ledges, 
                addNodes = list(
                  list(label = "Selected drug", shape = "diamond", color = "blue"),
                  list(label = "Drug", shape = "diamond", color = "skyblue"),
                  list(label = "Tumor", shape = "dot", color = "peachpuff")))
    
    
    return(network)
    
    
  }
}


#Plotly empty plot
ShinyDeepDR_EmptyPlot <- function(title = NULL){
  p <- plotly_empty(type = "scatter", mode = "markers") %>%
    config(
      displayModeBar = FALSE
    ) %>%
    layout(
      title = list(
        text = title,
        yref = "paper",
        y = 0.5
      )
    )
  return(p)
} 

#Created hyperlink in data table

ShinyDeepDR_CreateLink <- function(PubChem_CID) {
  paste0('<a href="https://pubchem.ncbi.nlm.nih.gov/compound/',PubChem_CID,'" target= "_blank">PubChem</a>')
}