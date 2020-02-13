
#' Creates a heatmap of the relative transmission probabilities
#'
#' The function \code{nbHeatmap} creates a heatmap of the transmission probabilities.
#' The rows are the possible infectors and the columns are the possible infectees both
#' ordered by \code{<dateVar>}. The darker the square the higher the probability that
#' the pair represented by that square is a transmission link. If a cluster method is specified 
#' using \code{clustMethod} and \code{cutoff}, then stars will be drawn in the squares of the
#' infectors in the top cluster.
#' 
#' Users have the option of specifying how the probabilities should be grouped into different
#' color shades through the argument \code{probBreaks}. The probabilities are split into groups by
#' using \code{probBreaks} as the \code{breaks} argument in \code{\link[base]{cut}} with the default options.
#' The length of the vector should be between 3 and 10 and the first element should be less than 0 and 
#' the last 1 so that all probabilities are guarenteed to be classified.
#' The colors are defined with the code \code{brewer.pal(length(probBreaks) - 1, "Blues")}
#' (where "Blues" is replaced by "Greys" if \code{blackAndWhite} is set to \code{TRUE}).
#' 
#' \strong{NOTE: This plot will take long to run and may not look good with
#'  larger outbreaks (>200 individuals)}
#' 
#' 
#' @param df The name of the dateset with transmission probabilities (column \code{pVar}),
#' individual IDs (columns \code{<indIDVar>.1} and \code{<indIDVar>.2}), and the dates of
#' observation (columns \code{<dateVar>.1} and \code{<dateVar>.2}).
#' @param indIDVar The name (in quotes) of the individual ID columns
#' (data frame \code{df} must have variables called \code{<indIDVar>.1}
#'  and \code{<indIDVar>.2}).
#' @param dateVar The name (in quotes) of the columns with the dates that the individuals are
#' observed (data frame \code{df} must have variables called \code{<dateVar>.1} and
#' \code{<dateVar>.2}).
#' @param pVar The name (in quotes) of the column with transmission probabilities.
#' @param clustMethod The method used to cluster the infectors; one of 
#' \code{"none", "n", "kd", "hc_absolute", "hc_relative"} where \code{"none"} or
#' not specifying a value means use all pairs with no clustering
#' (see \code{\link{clusterInfectors}} for detials on clustering methods).
#' @param cutoff The cutoff for clustering (see \code{\link{clusterInfectors}}).
#' @param blackAndWhite A logical. If \code{TRUE}, then the squares are colored in greyscale,
#' if \code{FALSE}, then the squares are colored with shades of blue.
#' @param probBreaks A numeric vector containing between 3 and 10 elements specifying the
#' boundaries used to classify the probabilities and color the squares.
#' The first element should be less than 0 and the last should be 1.
#' 
#' 
#' @examples
#' 
#' ## Heatmap with no clustering in color with the default probability breaks
#' par(mar = c(0, 0, 1, 0))
#' nbHeatmap(nbResults, indIDVar = "individualID", dateVar = "infectionDate",
#' pVar = "pScaled", clustMethod = "none")
#'
#' ## Adding stars for the top cluster, in black and white, changing the probability breaks
#' par(mar = c(0, 0, 1, 0))
#' nbHeatmap(nbResults, indIDVar = "individualID", dateVar = "infectionDate",
#'           pVar = "pScaled", clustMethod = "hc_absolute", cutoff = 0.05,
#'           blackAndWhite = TRUE, probBreaks = c(-0.01, 0.01, 0.1, 0.25, 0.5, 1))
#' 
#' 
#' @seealso \code{\link{nbProbabilities}} \code{\link{clusterInfectors}}
#' 
#' @export

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
  net.adj <- igraph::get.adjacency(net, attr = pVar, sparse = FALSE)
  
  #Marking cluster 1 with a * if clustering was specified
  net.trans <- igraph::get.adjacency(net, attr = "cluster", sparse = FALSE)
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
  
  pheatmap::pheatmap(t(net.adj), cluster_rows = FALSE, cluster_cols = FALSE,
           border_color = NA, show_rownames = FALSE, show_colnames = FALSE,
           col=RColorBrewer::brewer.pal(length(probBreaks) - 1, edgeCol), 
           breaks = probBreaks, display_numbers = t(net.trans), number_color = "white",
           fontsize_number = 7)
  
}




#' Creates a network of the relative transmission probabilities
#'
#' The function \code{nNetwork} creates a network of the transmission probabilities.
#' The nodes are the individuals and the edges represent possible transmission pairs.
#' The darker the edge, the higher the probability that the pair is a transmission link.
#' If a cluster method is specified using \code{clustMethod} and \code{cutoff}, only edges
#' that are in the top cluster of infectors will be drawn.
#' 
#' Users have the option of specifying how the probabilities should be grouped into different
#' color shades through the argument \code{probBreaks}. The probabilities are split into groups by
#' using \code{probBreaks} as the \code{breaks} argument in \code{\link[base]{cut}} with the default options.
#' The length of the vector should be between 3 and 10 and the first element should be less than 0 and 
#' the last 1 so that all probabilities are guarenteed to be classified.
#' The colors are defined with the code \code{brewer.pal(length(probBreaks) - 1, "Blues")}
#' (where "Blues" is replaced by "Greys" if \code{blackAndWhite} is set to \code{TRUE}).
#' 
#' 
#' @param df The name of the dateset with transmission probabilities (column \code{pVar}),
#' individual IDs (columns \code{<indIDVar>.1} and \code{<indIDVar>.2}), and the dates of
#' observation (columns \code{<dateVar>.1} and \code{<dateVar>.2}).
#' @param indIDVar The name (in quotes) of the individual ID columns
#' (data frame \code{df} must have variables called \code{<indIDVar>.1}
#'  and \code{<indIDVar>.2}).
#' @param dateVar The name (in quotes) of the columns with the dates that the individuals are
#' observed (data frame \code{df} must have variables called \code{<dateVar>.1} and
#' \code{<dateVar>.2}).
#' @param pVar The name (in quotes) of the column with transmission probabilities.
#' @param clustMethod The method used to cluster the infectors; one of 
#' \code{"none", "n", "kd", "hc_absolute", "hc_relative"} where \code{"none"} or
#' not specifying a value means use all pairs with no clustering
#' (see \code{\link{clusterInfectors}} for detials on clustering methods).
#' @param cutoff The cutoff for clustering (see \code{\link{clusterInfectors}}).
#' @param blackAndWhite A logical. If \code{TRUE}, then the edges are colored in greyscale,
#' if \code{FALSE}, then the edges are colored with shades of blue.
#' @param probBreaks A numeric vector containing between 3 and 10 elements specifying the
#' boundaries used to classify the probabilities and color the edges.
#' The first element should be less than 0 and the last should be 1.
#' 
#' 
#' @examples
#' 
#' ## Network of all pairs in color with the default probability breaks
#' par(mar = c(0, 0, 0.2, 0))
#' nbNetwork(nbResults, indIDVar = "individualID", dateVar = "infectionDate",
#' pVar = "pScaled", clustMethod = "none")
#'
#'## Adding stars for the top cluster, in black and white, changing the probability breaks
#' par(mar = c(0, 0, 0.2, 0))
#' nbNetwork(nbResults, indIDVar = "individualID", dateVar = "infectionDate",
#'           pVar = "pScaled", clustMethod = "hc_absolute", cutoff = 0.05,
#'           blackAndWhite = TRUE, probBreaks = c(-0.01, 0.01, 0.1, 0.25, 0.5, 1))
#' 
#' 
#' @seealso \code{\link{nbProbabilities}} \code{\link{clusterInfectors}}
#' 
#' @export


nbNetwork <- function(df, indIDVar, dateVar, pVar,
                      clustMethod = c("none", "n", "kd",
                                      "hc_absolute", "hc_relative"),
                      cutoff = NA, blackAndWhite = FALSE,
                      probBreaks = c(-0.01, 0.001, 0.005, 0.01,
                                     0.05, 0.1, 0.25, 0.5, 0.75, 1)){
  
  #Create a network of the probabilities
  net <- createNetwork(df, indIDVar = indIDVar, dateVar = dateVar, pVar  = pVar,
                      clustMethod = clustMethod, cutoff = cutoff,
                      probBreaks = probBreaks)
  
  #Network of top cluster
  net_top <- igraph::delete.edges(net, igraph::E(net)[igraph::E(net)$cluster == 2])
  table(igraph::E(net_top)$cluster)
  
  #Setting the edge colors to blue unless blackAndWhite is TRUE
  if(blackAndWhite == TRUE){
    edgeCol = "Greys"
  }else{
    edgeCol = "Blues"
  }

  graphics::plot(net_top, vertex.size = 7, vertex.label.cex = 0.7,
       vertex.color = "gray", vertex.frame.color = "dark gray",
       edge.width = 2, edge.arrow.size = 0.4,
       edge.color = RColorBrewer::brewer.pal(length(probBreaks) - 1,
                               edgeCol)[igraph::E(net_top)$pGroup])
  
}



## Function that creates a network from probabilities; called by nbHeatmap and nbNetwork
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
  
  
  #Checking that the named variables are in the data frame
  if(!indIDVar1 %in% names(df)){
    stop(paste0(indIDVar1, " is not in the data frame."))
  }
  if(!indIDVar2 %in% names(df)){
    stop(paste0(indIDVar2, " is not in the data frame."))
  }
  if(!dateVar1 %in% names(df)){
    stop(paste0(dateVar1, " is not in the data frame."))
  }
  if(!dateVar2 %in% names(df)){
    stop(paste0(dateVar2, " is not in the data frame."))
  }
  if(!pVar %in% names(df)){
    stop(paste0(pVar, " is not in the data frame."))
  }
  #Renaming probability to pScaled because I need to call it using $ notation
  df$pScaled <- df[, pVar]
  
  #Checking that the date variables are in a date form
  if(lubridate::is.Date(df[, dateVar1]) == FALSE &
     lubridate::is.POSIXt(df[, dateVar1]) == FALSE){
    stop(paste0(dateVar1, " must be either a date or a date-time (POSIXt) object."))
  }
  if(lubridate::is.Date(df[, dateVar2]) == FALSE &
     lubridate::is.POSIXt(df[, dateVar2]) == FALSE){
    stop(paste0(dateVar2, " must be either a date or a date-time (POSIXt) object."))
  }
  
  
  #If clustMethod is not specified, setting it to "none"
  if(length(clustMethod) > 1){
    clustMethod <- "none"
    warning("No clustMethod was provided so it was set to 'none'")
  }
  if(!clustMethod %in% c("none", "n", "kd", "hc_absolute", "hc_relative")){
    stop("clustMethod must be one of: none, n, kd, hc_absolute, hc_relative")
  }
  
  #Make sure cutoff is provided if clustMethod is not "none"
  if(clustMethod != "none" & is.null(cutoff)){
    stop("Please provide one or more cutoff values")
  }
  
  #Making sure probBreaks has the right form
  testValue <- probBreaks <= 1
  if(sum(testValue) != length(probBreaks)){
   stop("All values of probBreaks should be less than 1") 
  }
  if(probBreaks[1] > 0){
    probBreaks <- c(-0.01, probBreaks)
    message("First element of probBreaks is not negative so -0.01 was added to the beginning")
  }
  if(probBreaks[length(probBreaks)] != 1){
    probBreaks <- c(probBreaks, 1)
    message("Last element of probBreaks is not 1 so 1 was added to the end")
  }
  if(length(probBreaks) < 3 | length(probBreaks) > 10){
    stop("Please make sure probBreaks has between 3 and 10 elements")
  }
  
  
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
  net <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = T)
  #Adding a probability group based on probBreaks
  igraph::E(net)$pGroup <- cut(igraph::E(net)$pScaled, breaks = probBreaks,
                       labels = 1:(length(probBreaks) - 1))
  
  return(net)
}




#' Creates a plot of the effective reproductive number
#'
#' The function \code{plotRt} creates a plot of the effective reproductive number (Rt) over
#' the course of the outbreak. Using various options, the plot can include the overall average
#' Rt value for the outbreak and the confidence intervals.
#' 
#' The main input \code{rData} should be the output of \code{\link{estimateRt}} with the
#' time-level reproductive numbers, overall average, range used to calculate that average,
#' and time frame.
#' 
#' The options \code{includeRtCI} and \code{includeRtAvgCI} add confidence interval bounds
#' to the plot. If set to true, \code{rData} should be from a call of \code{\link{estimateRt}}
#' with \code{bootSamples > 0} so that confidence intervals are available.
#' If \code{includeRtAvgCI} is set to \code{TRUE}, a line for the point estimate of the average
#' Rt value will be drawn even if \code{includeRtAvg} is set to \code{FALSE}.
#' 
#' 
#' @param rData A list that is the output of \code{\link{estimateR}}. It should contain
#' the dataframes \code{RtDf}, \code{RtAvgDf}, and vectors \code{timeFrame} and \code{rangeForAvg}
#' @param includeRtAvg A logical. If TRUE, a horizontal line will be drawn for the average
#' Rt value over \code{rangeForAvg} and verticle lines will be drawn at the 
#' \code{rangeForAvg} values.
#' @param includeRtCI A logical. If TRUE, error bars will be added to the Rt values
#' representing the bootstrap confidence intervals.
#' @param includeRtAvgCI A logical. If TRUE, horizontal lines will be drawn around the Rt average
#' line representing the bootstrap confidence interval.
#' 
#' 
#' @examples
#' 
#' ## Use the nbResults data frame included in the package which has the results
#' # of the nbProbabilities() function on a TB-like outbreak.
#' 
#' ## Getting initial estimates of the reproductive number
#' # (ithout specifying nbResults and without confidence intervals)
#' rInitial <- estimateR(nbResults, dateVar = "infectionDate",
#'                indIDVar = "individualID", pVar = "pScaled",
#'                timeFrame = "months")
#'                
#' ## Finding the stable portion of the outbreak for rangeForAvg using the plot
#' plotRt(rInitial)
#' cut1 <- 25
#' cut2 <- 125
#' 
#' ## Finding the final reproductive number estimates with confidence intervals
#' # NOTE should run with bootSamples > 10.
#' rFinal <- estimateR(nbResults, dateVar = "infectionDate",
#'              indIDVar = "individualID", pVar = "pScaled",
#'              timeFrame = "months", rangeForAvg = c(cut1, cut2),
#'              bootSamples = 10, alpha = 0.05)
#' 
#' ## Ploting the final result              
#' plotRt(rFinal, includeRtAvg = TRUE, includeRtCI = TRUE, includeRtAvgCI = TRUE)
#' 
#' @seealso \code{\link{nbProbabilities}} \code{\link{estimateR}}
#' 
#' @export


plotRt <- function(rData, includeRtAvg = FALSE,
                   includeRtCI = FALSE, includeRtAvgCI = FALSE){
  
  #Checking to make sure rData has the right form
  if(class(rData) != "list"){
    stop("The rData argument should be the list output from the function estimateR")
  }
  if(identical(names(rData), c("RiDf", "RtDf", "RtAvgDf", "timeFrame", "rangeForAvg")) == FALSE){
    stop("The rData argument should be the list output from the function estimateR")
  }
  
  #Translating the timeFrame into labels for the axes
  xlabel <- gsub("s$", "", Hmisc::capitalize(rData$timeFrame))
  if(rData$timeFrame == "days"){
    ylabel <- "Daily"
  }else{
    ylabel <- paste0(xlabel, "ly")
  }
  
  #Base plot of reproductive number over time
  p <- ggplot2::ggplot(data = rData$RtDf, ggplot2::aes(x = rData$RtDf$timeRank, y = rData$RtDf$Rt)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(name = paste0(ylabel, " Effective Reproductive Number")) + 
    ggplot2::scale_x_continuous(name = paste0(xlabel, " of Observation")) +
    ggplot2::theme_bw() 
  
  #Adding confidence interval for Rt
  if(includeRtCI == TRUE){
    
    if(!"ciLower" %in% names(rData$RtDf)){
      stop("Please provide a rData list that has confidence intervals")
    }
    p <- p +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = rData$RtDf$ciLower, ymax = rData$RtDf$ciUpper),
                    width = 0.1, color = "darkgrey")
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
      ggplot2::geom_vline(ggplot2::aes(xintercept = cut1), linetype = 3, size = 0.7) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = cut2), linetype = 3, size = 0.7) +
      ggplot2::geom_hline(data = rData$RtAvgDf, ggplot2::aes(yintercept = rData$RtAvg$RtAvg),
                          size = 0.7)
    
    #Adding confindeince interval for RtAvg
    if(includeRtAvgCI == TRUE){
      
      if(!"ciLower" %in% names(rData$RtAvgDf)){
        stop("Please provide a rData list that has confidence intervals")
      }
      
      p <- p +
        ggplot2::geom_hline(data = rData$RtAvgDf, ggplot2::aes(yintercept = rData$RtAvgDf$ciLower),
                   linetype = 2, size = 0.7, color = "darkgrey") +
        ggplot2::geom_hline(data = rData$RtAvgDf, ggplot2::aes(yintercept = rData$RtAvgDf$ciUpper),
                   linetype = 2, size = 0.7, color = "darkgrey") 
    }
  }
  
  print(p)
}


