


nbHeatmap <- function(df, indIDVar, dateVar, pVar,
                      clustMethod = c("none", "n", "kd",
                                      "hc_absolute", "hc_relative"),
                      cutoff = NA, blackAndWhite = FALSE,
                      probBreaks = c(-0.01, 0.001, 0.005, 0.01,
                                     0.05, 0.1, 0.25, 0.5, 0.75, 1)){
  
  #Create a network of the probabilities
  net <- createNetwork(df, indIDVar = indIDVar, dateVar = dateVar, pVar  = pVar,
                       clustMethod = clustMethod, cutoff = cutoff, probBreaks = probBreaks)
  
  #First get adjacency version of the network using pScaled
  net.adj <- get.adjacency(net, attr = pVar, sparse = FALSE)
  
  #Marking cluster 1 with a * if clustering was specified
  net.trans <- get.adjacency(net, attr = "cluster", sparse = FALSE)
  if(clustMethod != "none"){
    net.trans <- ifelse(net.trans == 1, "*", "") 
  }else{
    net.trans <- ifelse(net.trans == 1, "", "") 
  }
  
  #Setting the edge colors to blue unless blackAndWhite is TRUE
  if(blackAndWhite == TRUE){
    edgeCol = "Greys"
  }else{
    edgeCol = "Blues"
  }
  
  par(mar = c(0, 0, 1, 0))
  pheatmap(t(net.adj), cluster_rows = FALSE, cluster_cols = FALSE,
           border_color = NA, show_rownames = FALSE, show_colnames = FALSE,
           col=brewer.pal(length(probBreaks) - 1, edgeCol), breaks = probBreaks,
           display_numbers = t(net.trans), number_color = "white",
           fontsize_number = 7)
  
}


nbNetwork <- function(df, indIDVar, dateVar, pVar,
                      clustMethod = c("none", "n", "kd",
                                      "hc_absolute", "hc_relative"),
                      cutoff = NA, blackAndWhite = FALSE,
                      probBreaks = c(-0.01, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)){
  
  #Create a network of the probabilities
  net <- createNetwork(df, indIDVar = indIDVar, dateVar = dateVar, pVar  = pVar,
                      clustMethod = clustMethod, cutoff = cutoff, probBreaks = probBreaks)
  
  #Network of top cluster
  net_top <- delete.edges(net, E(net)[cluster == 2])
  
  #Setting the edge colors to blue unless blackAndWhite is TRUE
  if(blackAndWhite == TRUE){
    edgeCol = "Greys"
  }else{
    edgeCol = "Blues"
  }

  par(mar = c(0, 0, 0.2, 0))
  plot(net_top, vertex.size = 7, vertex.label.cex = 0.7,
       vertex.color = "gray", vertex.frame.color = "dark gray",
       edge.width = 2, edge.arrow.size = 0.4,
       edge.color = brewer.pal(length(probBreaks) - 1,
                               edgeCol)[E(net_top)$pGroup])
  
}


createNetwork <- function(df, indIDVar, dateVar, pVar,
                            clustMethod = c("none", "n", "kd",
                                            "hc_absolute", "hc_relative"),
                            cutoff = NA, probBreaks = c(-0.01, 0.001, 0.005, 0.01,
                                                        0.05, 0.1, 0.25, 0.5, 0.75, 1)){
  
  df <- as.data.frame(df)
  #Creating variables with the individual ID and date variables
  indIDVar1 <- paste0(indIDVar, ".1")
  indIDVar2 <- paste0(indIDVar, ".2")
  dateVar1 <- paste0(dateVar, ".1")
  dateVar2 <-  paste0(dateVar, ".2")
  #Renaming probability to pScaled because need to call it using $ notation
  df$pScaled <- df[, pVar]
  
  #If clustMethod is "none" setting assigning all pairs to cluster 1 or clustering
  #the probabilities based on clustMethod and cutoff, then restricting to
  #just the top cluster to be used for estimation
  if(clustMethod == "none"){
    clustRes <- df
    clustRes$cluster <- 1
  }else{
    clustRes <- clusterInfectors(df, indIDVar = indIDVar, pVar = pVar,
                                 clustMethod = clustMethod, cutoff = cutoff)
  }
  topClust <- clustRes[clustRes$cluster == 1, ]
  
  
  ## Setting up the network ##

  #Finding the individual ids and dates for each individual
  indIDs1 <- df[!duplicated(df[, indIDVar1]),
                c(indIDVar1, dateVar1)]
  names(indIDs1) <- c(indIDVar, dateVar)
  indIDs2 <- df[!duplicated(df[, indIDVar2]),
                c(indIDVar2, dateVar2)]
  names(indIDs2) <- c(indIDVar, dateVar)
  
  #Creating a dataframe of nodes
  nodes <- rbind(indIDs1, indIDs2)
  nodes <- nodes[!duplicated(nodes[, indIDVar]), ]
  nodes <- nodes[order(nodes[, dateVar]), ]
  nodes[, dateVar] <- as.character(dateVar)
  
  #Creating a dataframe of edges
  #Reording the columns so individualIDs are first
  #Also removing label column because it causes issues
  notIDNames <- names(clustRes)[!names(clustRes) %in% c(indIDVar1, indIDVar2, "label")]
  namesReordered <- c(indIDVar1, indIDVar2, notIDNames)
  edges <- clustRes[order(clustRes[, pVar]), namesReordered]
  
  #Setting the date variables to be characters to avoid warning
  for(i in 1:ncol(edges)){
    if("POSIXct" %in% class(edges[, i])){
      edges[, i] <- as.character(edges[, i])
    }
  }
  
  #Creating the network
  net <- graph_from_data_frame(d = edges, vertices = nodes, directed = T)
  #Adding a probability group based on probBreaks
  E(net)$pGroup <- cut(E(net)$pScaled, breaks = probBreaks,
                       labels = 1:(length(probBreaks) - 1))
  
  return(net)
}



rtPlot <- function(df){
  
}