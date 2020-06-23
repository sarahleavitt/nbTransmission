
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
#' greater than the lowest 0 region are assigned to the high probability cluster. If the density of the
#' probabilities does not drop to 0 then all infectors are assigned to the low probability cluster 
#' (indicating no real clustering).
#' 
#' If \code{clustMethod == "hc_absolute"} or \code{clustMethod == "hc_relative"}, then
#' hierarchical clustering with minimum distance is used to split the possible infectors
#' into two clusters. This method functionally splits the infectors by the largest gap
#' in their probabilities.
#' 
#' Then if \code{clustMethod == "hc_absolute"}, those infectees
#' where the gap between the two clusters is less than \code{cutoff} have all of their
#' possible infectors reassigned to the low probability cluster (indicating no real clustering).
#' If \code{clustMethod == "hc_relative"}, then all infectees where the gap between the two
#' clusters is less than \code{cutoff} times the second largest gap in probabilities
#' are reassigned to the low probability cluster (indicating no real clustering).
#' 
#' 
#' @param df The name of the dateset with transmission probabilities (column \code{pVar}),
#' individual IDs (columns \code{<indIDVar>.1} and \code{<indIDVar>.2}).
#' @param indIDVar The name (in quotes) of the individual ID columns
#' (data frame \code{df} must have variables called \code{<indIDVar>.1}
#'  and \code{<indIDVar>.2}).
#' @param pVar The name (in quotes) of the column with transmission probabilities.
#' @param clustMethod The method used to cluster the infectors; 
#' one of \code{"n", "kd", "hc_absolute", "hc_relative"} (see details).
#' @param cutoff The cutoff for clustering (see details).
#' 
#'
#' @return The original data frame (\code{df}) with a new column called \code{cluster}
#' which is a factor variable with value \code{1} if the infector is in the high probability cluster
#' or \code{2} if the infector is in the low probability cluster.
#' 
#' 
#' @examples
#' 
#' ## Use the nbResults data frame included in the package which has the results
#' ## of the nbProbabilities() function on a TB-like outbreak.
#' 
#' ## Clustering using top n
#' # High probability cluster includes infectors with highest 3 probabilities
#' clust1 <- clusterInfectors(nbResults, indIDVar = "individualID", pVar = "pScaled",
#'                            clustMethod = "n", cutoff = 3)
#' table(clust1$cluster)
#' 
#' ## Clustering using hierarchical clustering
#'
#' # Cluster all infectees, do not force gap to be certain size
#' clust2 <- clusterInfectors(nbResults, indIDVar = "individualID", pVar = "pScaled",
#'                             clustMethod = "hc_absolute", cutoff = 0)
#' table(clust2$cluster)
#' 
#' ## NOT RUN ##
#' # # Absolute difference: gap between top and bottom clusters is more than 0.05
#' # clust3 <- clusterInfectors(nbResults, indIDVar = "individualID", pVar = "pScaled",
#' #                            clustMethod = "hc_absolute", cutoff = 0.05)
#' # table(clust3$cluster)
#'
#' ## NOT RUN ##
#' # # Relative difference: gap between top and bottom clusters is more than double any other gap
#' # clust4 <- clusterInfectors(nbResults, indIDVar = "individualID", pVar = "pScaled",
#' #                            clustMethod = "hc_relative", cutoff = 2)
#' # table(clust4$cluster)
#'
#' ## NOT RUN ##
#' # ## Clustering using kernel density estimation
#' # # Using a small binwidth of 0.01
#' # clust5 <- clusterInfectors(nbResults, indIDVar = "individualID", pVar = "pScaled",
#' #                            clustMethod = "kd", cutoff = 0.01)
#' # table(clust5$cluster)
#' 
#' @seealso \code{\link{nbProbabilities}}
#' 
#' @export



clusterInfectors <- function(df, indIDVar, pVar,
                             clustMethod = c("n", "kd", "hc_absolute", "hc_relative"),
                             cutoff){
  
  df <- as.data.frame(df)
  #Creating variables with the individual ID
  indIDVar1 <- paste0(indIDVar, ".1")
  indIDVar2 <- paste0(indIDVar, ".2")
  
  #Checking that the named variables are in the data frame
  if(!indIDVar1 %in% names(df)){
    stop(paste0(indIDVar1, " is not in the data frame."))
  }
  if(!indIDVar2 %in% names(df)){
    stop(paste0(indIDVar2, " is not in the data frame."))
  }
  if(!pVar %in% names(df)){
    stop(paste0(pVar, " is not in the data frame."))
  }
  
  #Making sure clustMethod is correctly specified
  if(length(clustMethod) > 1){
    stop("Please provide a clustering method")
  }
  else if(!clustMethod %in% c("n", "kd", "hc_absolute", "hc_relative")){
    stop("clustMethod must be one of: n, kd, hc_absolute, hc_relative")
  }
  
  
  #Ranking the probabilities for each possible infector
  #Ties are set to the minimum rank of that group
  df <- df[order(df[, indIDVar2], -df[, pVar]), ]
  df$pRank <- stats::ave(-df[, pVar], df[, indIDVar2], 
                             FUN = function(x){
                               rank(x, ties.method = "min") 
                             })
  
  #Adding the number of infectors
  nInf <- as.data.frame(table(df[, indIDVar2]))
  names(nInf) <- c(indIDVar2, "nInfectors")
  df2 <- merge(df, nInf, by = indIDVar2, all = TRUE)
  
  #Splitting the data frame to those who have one infector and multiple infectors
  multInf <- df2[df2$nInfectors > 1, ]
  oneInf <- df2[df2$nInfectors == 1, ]
  
  ## Using a constant number of infectors ##
  if(clustMethod == "n"){
    
    clustRes <- multInf
    clustRes$cluster <- ifelse(clustRes$pRank <= cutoff, 1, 2)
    
  }
  
  ## Using kernel density estimation ##
  if(clustMethod == "kd"){
    
    clustRes1 <- dplyr::group_by(multInf, !!rlang::sym(indIDVar2))
    clustRes <- dplyr::group_modify(clustRes1, ~ findClustersKD(.x, pVar = pVar,
                                                     cutoff = cutoff))

  }
  
  ## Using hierarchical clustering ##
  if(grepl("hc", clustMethod)){
    
    clustRes1 <- dplyr::group_by(multInf, !!rlang::sym(indIDVar2))
    clustRes <- dplyr::group_modify(clustRes1, ~ findClustersHC(.x, pVar,
                                                                 cutoff = cutoff,
                                                                 clustMethod = clustMethod))
  
  }
  
  #Removing tibble formatting
  clustRes <- as.data.frame(dplyr::ungroup(clustRes))
  
  #Combining the clustering for those with more than one infector with those
  #who have one infector
  if(nrow(oneInf) > 0){
    #If there is one infector, it should be in the top cluster
    oneInf$cluster <- 1
    clustRes2 <- rbind(oneInf, clustRes)
  }else{
    clustRes2 <- clustRes
  }
  
  #Making the clusters a factor variable
  clustRes2$cluster <- factor(clustRes2$cluster, levels = c(1, 2))
  #Removing nInfectors variable
  clustRes2$nInfectors <- NULL
  
  return(clustRes2)
}





## Function to find clusters using kernel density estimation ##

findClustersKD <- function(df, pVar, cutoff = 0.05, minGap = 0, plot = FALSE,
                           colors = c("#00BFC4", "#F8766D"), size){
  
  df <- as.data.frame(df)
  df <- df[order(df[, pVar]),]
  df$cluster <- NA
  
  tryCatch({
    
    #Estimating the density for the probabilities using the cutoff as the binwidth
    d <- stats::density(df[, pVar], bw = cutoff, kernel = "rectangular", from = 0)
    
    #Finding the indices where the density drops to 0 indicating a split
    mindf <- cbind.data.frame(index = which(d$y < 0.00001),
                              minx = d$x[which(d$y < 0.00001)])
    
    #Finding the difference in the indices for each time the density goes to 0
    mindf$xdiff = mindf$index - dplyr::lag(mindf$index)
    mindf$xdiff[is.na(mindf$xdiff)] <- 1
    #Restricting to the range of the original probabilities
    mindf <- mindf[mindf$minx > min(df[, pVar]) & 
                     mindf$minx < max(df[, pVar]), ]
    
    #Finding all of the separate regions of x where the density goes to 0
    #If xdiff > 1 then this indicates a sperate region of 0 density
    if(nrow(mindf) > 1){
      
      region <- 1
      for(i in 1:nrow(mindf)){
        row <- mindf[i, ]
        if(mindf[i, "xdiff"] == 1){
          region = region
        }else{
          region = region + 1
        }
        mindf[i, "region"] <- region
      }
      
      #Finding the length of these minimum regions
      groupMinsL <- by(mindf,
                   INDICES = list(mindf$region),
                   FUN = function(x){
                     data.frame("region" = unique(x$region),
                                "lower" = min(x$minx, na.rm = TRUE),
                                "upper" = max(x$minx, na.rm = TRUE))
                   })
      groupMins <- do.call(rbind, groupMinsL)
      #Making sure very small lengths are 0
      groupMins$length <- ifelse(groupMins$upper - groupMins$lower < 0.000001, 0,
                                 groupMins$upper - groupMins$lower)
      
      #Restricting to regions that are longer than the minGap
      groupMins <- groupMins[groupMins$length > minGap, ]
      
      #Finding the x value at the upper bound of the first region with length more
      #than minGap where the density goes to 0 and call it lowestMin.
      #If there is no such region set lowestMin to the highest probability.
      if(nrow(groupMins) != 0){
        lowestMin <- min(groupMins$upper)
      }else{
        lowestMin = max(df[, pVar])
      }
    #If there are no places that the density goes to 0 (mindf is empty)
    #Then set lowestMin to the highest probability
    }else{
      lowestMin = max(df[, pVar])
    }
    
    #Create the clustering value such that an infector is in cluster 1 if it
    #has a probability above the lowestMin (above the region where the density goes to 0)
    #and is in cluster 2 if the probability is below that region or if there is no
    #such region, all infectors are in cluster 2.
    df$cluster <- ifelse(df[, pVar] > lowestMin, 1, 2)
    
    
    if(plot == TRUE){
      densitydf <- cbind.data.frame(x = d$x, y = d$y)
      
      p <- ggplot2::ggplot(data = df) +
        ggplot2::geom_histogram(ggplot2::aes(x = df[, pVar], fill = factor(df$cluster, levels = c(1, 2))), bins = 20) +
        ggplot2::geom_line(data = densitydf, ggplot2::aes(x = densitydf$x, y = densitydf$y),
                           color = "black", alpha = 0.5) +
        ggplot2::xlab("Relative Probability") +
        ggplot2::ylab("Count") +
        ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
        ggplot2:: theme_bw(base_size = size) +
        ggplot2::theme(legend.position = "none")
      
      print(p)
      return(p)
    }
    
  }, error = function(e){
    print(nrow(df))
    cat("ERROR: ", conditionMessage(e), "\n")})
  
  return(df)
}





## Function to find clusters using hierarchical clustering ##

findClustersHC <- function(df, pVar, cutoff = 0.05,
                           clustMethod = c("hc_absolute", "hc_relative")){
  
  df <- as.data.frame(df[order(-df[, pVar]), ])
  
  #Clustering the infectors using hierarchical cluster with minimum distance
  hclustS <- stats::hclust(stats::dist(df[, pVar]), method = "single")
  hclustScut <- stats::cutree(hclustS, 2)
  df$cluster <- hclustScut
  
  #Finding the boundaries of the two clusters:
  #Minimum value in the top cluster
  minC1 <- min(df[df$cluster == 1, pVar])
  #Maximum value in the bottom cluster
  maxC2 <- max(df[df$cluster == 2, pVar])
  #Finding the gap between the two clusters
  gap <- minC1 - maxC2
  
  ## Alternatively ##
  #gap <- hclustS$height[length(hclustS$height)]
  
  #Determining which cases should have clusters based on the absolute
  #size of the gap between clusters
  if(clustMethod == "hc_absolute"){
    
    if(gap < cutoff){
      df$cluster <- 2
      }
  }
  
  #Determining which cases should have clusters based on the relative
  #size of the gap between clusters - gap between clusters needs to be
  #cutoff times greater then the next largest gap in probabilities
  if(clustMethod == "hc_relative"){

    allGaps <- diff(df[, pVar])
    secondGap <- abs(sort(allGaps)[2])
    
    if(is.na(secondGap) | gap < cutoff*secondGap){
      df$cluster <- 2
      }
  }
  return(df)
}

