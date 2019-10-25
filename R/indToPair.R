
#' Individual to Pair
#'
#' Transforms a dataset of individuals to a dataset of pairs
#' 
#' Add details section
#'
#' @param indData An individual-level dataframe.
#' @param indIDVar The variable name (in quotes) of the individual ID variable.
#' @param separator The character to be used to separate the individual IDs when creating
#' the pairID
#' @param dateVar The variable name (in quotes) of the dates that the individuals are observed
#' (optional unless ordered = TRUE). If supplied, the time difference in days and years between
#' individuals will be calculated.
#' @param ordered A logical indicating if a set of ordered pairs should be returned
#' (dateVar.1 before dateVar.2). If FALSE a dataframe of all pairs will be returned
#'
#' @return A dataframe of either all possible pairs of individuals (ordered = FALSE) or ordered
#' pairs of individuals (ordered = TRUE). The dataframe will have all of the original variables
#' with suffixes ".1" and ".2" corresponding to the original values of 
#' \code{indIDVar.1} and \code{indIDVar.2}.
#' 
#' Added to the dataframe will be a column called \code{pairID} which is \code{indIDVar.1}
#' and \code{indIDVar.2} separated by \code{separator}.
#' 
#' If dateVar is provided the dataframe will also include variables \code{timeDiff} and 
#' \code{timeDiffY} giving the difference of time in days and years respectively between the value
#' of \code{dateVar} for \code{indIDVar.1} and \code{indIDVar.2}
#'
#' @examples
#' #Insert example here
#' 
#' @export
#' 



indToPair <- function(indData, indIDVar, separator = "_", dateVar = NULL, ordered = FALSE){
  
  indData <- as.data.frame(indData)
  #Finding all pairs of IDs (order matters)
  pairs <- expand.grid(indData[, indIDVar], indData[, indIDVar])
  
  #Removing pairs of the same individual
  pairs2 <- pairs[pairs$Var1 != pairs$Var2, ]
  
  #Creating an pairID that combines the individualIDs with an underscore
  pairs2$pairID <- paste(pairs2$Var1, pairs2$Var2, sep = separator)
  
  #Renaming the dataset with the individual ID
  names(pairs2) <- c(paste0(indIDVar, ".1"), paste0(indIDVar, ".2"), "pairID")
  
  #Merging back with the rest of the variables
  pairData1 <- merge(pairs2, indData, by.x = paste0(indIDVar, ".1"), by.y = indIDVar,
                     all = TRUE)
  pairData2 <- merge(pairData1, indData, by.x = paste0(indIDVar, ".2"), by.y = indIDVar,
                     all = TRUE, suffixes = c(".1", ".2"))
  
  if(!is.null(dateVar)){
    pairData2$timeDiff <- as.numeric(difftime(pairData2[, paste0(dateVar, ".2")],
                                              pairData2[, paste0(dateVar, ".1")],
                                              units = "days"))
    pairData2$timeDiffY <- pairData2$timeDiff / 365
  }
  
  if(ordered == TRUE & is.null(dateVar)){
    stop("If ordered = TRUE, then dateVar must be provided")
  }
  
  if(ordered == TRUE){
    orderedData <- pairData2[pairData2$timeDiff > 0, ]
    return(orderedData)
  }else{
    return(pairData2)
  }
}


