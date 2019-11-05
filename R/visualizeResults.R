
#' Clusters the infectors based on their transmission probabilities
#'
#' The function \code{clusterInfectors} uses either kernel density estimation or
#' hierarchical clustering to cluster the infectors for each infectee. This clustering
#' provides a way to separate out the few top possible infectors for each infectee
#' if there is such a cluster.
#' 
#' This function provides a way to find the most likely infectors for each infectee
#' using various clustering methods indicated by the \code{clustmethod}.
#' The methods can be one of \code{c("n", "kd", "hc_constant", "hc_relative")}.
#' 
#' If \code{clustMethod == "n"} then this function simply assigns the top n possible 
#' infectors in the top cluster where n is defined by the value of \code{cutoff}.
#' 
#' If \code{clustMethod == "kd"} then kernel density estimation is used to split the infectors.
#' The density for the probabilities for all infectors is estimated using a binwidth defined
#' by the value of \code{cutoff}. If the density is made up of at least two separate curves
#' (separated by a region where the density drops to 0) then the infectors with probabilities
#' greater than the lowest 0 region are assigned into the top cluster. If the density of the
#' probabilities does not drop to 0 then all infectors are assigned into the bottom cluster.
#' 
#' If \code{clustMethod == "hc_absolute"} or \code{clustMethod == "hc_relative"}, then
#' hierarchical clustering with minimum distance is used to split the possible infectors
#' into two clusters. This method functionally splits the infectors by the largest gap
#' in their probabilities.
#' 
#' Then if \code{clustMethod == "hc_absolute"}, those infectees
#' where the gap between the two clusters is less than \code{cutoff} have all of their
#' possible infectors reassigned to the bottom cluster (indicating no real clustering).
#' If \code{clustMethod == "hc_relative"}, then all infectees where the gap between the two
#' clusters is less than \code{cutoff} times the second largest gap in probabilities
#' are reassigned to the bottom cluster (indicating no real clustering).
#' 
#' 
#' @param df The name of the dateset with transmission probabilities (column \code{pVar}),
#' individual IDs (columns \code{<indIDVar>.1} and \code{<indIDVar>.2}).
#' @param indIDVar The name (in quotes) of the individual ID columns
#' (data frame \code{df} must have variables called \code{<indIDVar>.1}
#'  and \code{<indIDVar>.2}).
#' @param pVar The name (in quotes) of the column with transmission probabilities.
#' @param clustMethod The method used to cluster the infectors (see details).
#' @param cutoff The cutoff for clustering (see details).
#' 
#'
#' @return The original data frame (\code{df}) with a new column called \code{cluster}
#' which is a factor variable with value \code{1} if the infector is in the top cluster
#' or \code{2} if the infector is in the bottom cluster.
#' 
#' 
#' @examples
#' 
#' ## Use the pairData dataset which represents a TB-like outbreak
#' # First create a dataset of ordered pairs
#' orderedPair <- pairData[pairData$infectionDiffY > 0, ]
#' 
#' ## Create a variable called snpClose that will define probable links
#' # (<3 SNPs) and nonlinks (>12 SNPs) all pairs with between 2-12 SNPs
#' # will not be used to train.
#' orderedPair$snpClose <- ifelse(orderedPair$snpDist < 3, TRUE,
#'                         ifelse(orderedPair$snpDist > 12, FALSE, NA))
#' table(orderedPair$snpClose)
#' 
#' ## Running the algorithm
#' # NOTE should run with nReps > 1
#' covariates = c("Z1", "Z2", "Z3", "Z4", "timeCat")
#' resGen <- nbProbabilities(orderedPair = orderedPair,
#'                             indIDVar = "individualID",
#'                             pairIDVar = "pairID",
#'                             goldStdVar = "snpClose",
#'                             covariates = covariates,
#'                             label = "SNPs", l = 1,
#'                             n = 10, m = 1, nReps = 1)
#'                             
#' ## Merging the probabilities back with the pair-level data
#' allProbs <- merge(resGen[[1]], orderedPair, by = "pairID", all = TRUE)
#' 
#' ## Clustering using top n
#' # Top cluster includes infectors with highest 3 probabilities
#' clust1 <- clusterInfectors(allProbs, indIDVar = "individualID", pVar = "pScaled",
#'                            clustMethod = "n", cutoff = 3)
#' table(clust1$cluster)
#' 
#' ## Clustering using hierarchical clustering
#'
#' # Cluster all infectees, do not force gap to be certain size
#' clust2 <- clusterInfectors(allProbs, indIDVar = "individualID", pVar = "pScaled",
#'                            clustMethod = "hc_absolute", cutoff = 0)
#' table(clust2$cluster)
#' 
#' # Absolute difference: gap between top and bottom clusters is more than 0.05
#' clust3 <- clusterInfectors(allProbs, indIDVar = "individualID", pVar = "pScaled",
#'                            clustMethod = "hc_absolute", cutoff = 0.05)
#' table(clust3$cluster)
#'
#' # Relative difference: gap between top and bottom clusters is more than double any other gap
#' clust4 <- clusterInfectors(allProbs, indIDVar = "individualID", pVar = "pScaled",
#'                            clustMethod = "hc_relative", cutoff = 2)
#' table(clust4$cluster)
#'
#' ## Clustering using kernel density estimation
#' # Using a small binwidth of 0.01
#' clust5 <- clusterInfectors(allProbs, indIDVar = "individualID", pVar = "pScaled",
#'                            clustMethod = "kd", cutoff = 0.01)
#' table(clust5$cluster)
#' 
#' @seealso \code{\link{nbProbabilities}}
#' 
#' @export

nbHeatmap <- function(df, indIDVar, dateVar, pVar,
                      clustMethod = c("none", "n", "kd",
                                      "hc_absolute", "hc_relative"),
                      cutoff = NA, blackAndWhite = FALSE,
                      probBreaks = c(-0.01, 0.001, 0.005, 0.01,
                                     0.05, 0.1, 0.25, 0.5, 0.75, 1)){
  
  #If clustMethod is not specified, setting it to "none"
  if(length(clustMethod) > 1){
    clustMethod <- "none"
  }
  
  #Create a network of the probabilities
  net <- createNetwork(df, indIDVar = indIDVar, dateVar = dateVar, pVar  = pVar,
                       clustMethod = clustMethod, cutoff = cutoff, probBreaks = probBreaks)
  
  #First get adjacency version of the network using pScaled
  net.adj <- get.adjacency(net, attr = pVar, sparse = FALSE)
  
  #Marking cluster 1 with a * if clustering was specified
  net.trans <- get.adjacency(net, attr = "cluster", sparse = FALSE)
  if(clustMethod == "none"){
    net.trans <- ifelse(net.trans == 1, "", "") 
  }else{
    net.trans <- ifelse(net.trans == 1, "*", "") 
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
  
  #Creating a data frame of nodes
  nodes <- rbind(indIDs1, indIDs2)
  nodes <- nodes[!duplicated(nodes[, indIDVar]), ]
  nodes <- nodes[order(nodes[, dateVar]), ]
  nodes[, dateVar] <- as.character(dateVar)
  
  #Creating a data frame of edges
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



rtPlot <- function(rData, includeRtAvg = TRUE,
                   includeRtCI = TRUE, includeRtAvgCI = TRUE){
  
  #Translating the timeFrame into labels for the axes
  xlabel <- gsub("s$", "", Hmisc::capitalize(rData$timeFrame))
  if(timeFrame == "days"){
    ylabel <- "Daily"
  }else{
    ylabel <- paste0(xlabel, "ly")
  }
  
  #Base plot of reproductive number over time
  p <- ggplot(data = rData$RtDf, aes(x = rData$RtDf$timeRank, y = rData$RtDf$Rt)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(name = paste0(ylabel, " Effective Reproductive Number")) + 
    scale_x_continuous(name = paste0(xlabel, " of Observation")) +
    theme_bw()
  
  #Adding confidence interval for Rt
  if(includeRtCI == TRUE){
    
    if(!"ciLower" %in% names(rData$RtDf)){
      stop("Please provide a rData list that has confidence intervals")
    }
    p <- p +
      geom_errorbar(aes(ymin = rData$RtDf$ciLower, ymax = rData$RtDf$ciUpper),
                    width = 0.1, color = "dark grey")
  }
  
  #Adding horizontal line for Rt and vertical lines showing range included in average
  if(includeRtAvg == TRUE | includeRtAvgCI == TRUE){
    
    if(is.null(rData$rangeForAvg)){
      cut1 <- 0
      cut2 <- max(rData$RtDf$timeRank)
    }else{
      cut1 <- rData$rangeForAvg[1]
      cut2 <- rData$rangeForAvg[2]
    }
    
    p <- p +    
      geom_vline(aes(xintercept = cut1), linetype = 3, size = 0.7) +
      geom_vline(aes(xintercept = cut2), linetype = 3, size = 0.7) +
      geom_hline(data = rFinal[[3]], aes(yintercept = RtAvg), size = 0.7)
    
    #Adding confindeince interval for RtAvg
    if(includeRtAvgCI == TRUE){
      
      if(!"ciLower" %in% names(rData$RtAvgDf)){
        stop("Please provide a rData list that has confidence intervals")
      }
      
      p <- p +
        geom_hline(data = rData$RtAvgDf, aes(yintercept = rData$RtAvgDf$ciLower),
                   linetype = 2, size = 0.7, color = "dark grey") +
        geom_hline(data = rData$RtAvgDf, aes(yintercept = rData$RtAvgDf$ciUpper),
                   linetype = 2, size = 0.7, color = "dark grey") 
    }
  }
  
  print(p)
}


