
#' Individual to Pair
#'
#' Transforms a dataset of individuals to a dataset of pairs
#' 
#' The function \code{indToPair} takes a dataset of observations (such as individuals in an infectious
#' disease outbreak) and transforms it into a dataset of pairs. The function requires an id column:
#' \code{indIDVar} to identify the individual observations. The resulting pair-level dataframe will
#' have a \code{pairID} column which combines the individual IDs for that pair.
#' 
#' The function can either output all possible pairs (\code{ordered = FALSE}) or only ordered pairs 
#' (\code{ordered = TRUE}) where the ordered is determined by a date variable (\code{dateVar}).
#' If \code{ordered = FALSE}, then the \code{dateVar} argument is option and if provided,
#' the output will include the time difference in days and years between the individuals in the pair.
#' If \code{orded = TRUE}, then \code{dateVar} must be provide.
#'
#' @param indData An individual-level dataframe.
#' @param indIDVar The name (in quotes) of the column with the individual ID.
#' @param separator The character to be used to separate the individual IDs when creating
#' the pairID.
#' @param dateVar The name (in quotes) of the column with the dates that the individuals are observed
#' (optional unless ordered = TRUE). This column must be a date or date-time (POSIXt) object.
#' If supplied, the time difference in days and years between individuals will be calculated.
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
#' If dateVar is provided the dataframe will also include variables \code{<dateVar>.Diff} and 
#' \code{<dateVar>.DiffY} giving the difference of time in days and years respectively between the value
#' of \code{dateVar} for \code{indIDVar.1} and \code{indIDVar.2}
#'
#' @examples
#' ## Create a dataset of all pairs with no date variable
#' pairU <- indToPair(indData = indData, indIDVar = "individualID")
#' 
#' ## Create a dataset of all pairs with a date variable
#' pairUD <- indToPair(indData = indData, indIDVar = "individualID",
#'                       dateVar = "infectionDate")
#' 
#' ## Create a dataset of ordered pairs
#' pairO <- indToPair(indData = indData, indIDVar = "individualID",
#'                      dateVar = "infectionDate", ordered = TRUE)
#' 
#' @export
#' 



indToPair <- function(indData, indIDVar, separator = "_", dateVar = NULL, ordered = FALSE){
  
  indData <- as.data.frame(indData)
  
  #Checking that the named variables are in the data frame
  if(!indIDVar %in% names(indData)){
    stop(paste0(indIDVar, " is not in the data frame."))
  }
  if(!is.null(dateVar)){
    if(!dateVar %in% names(indData)){
      stop(paste0(dateVar, " is not in the data frame.")) 
    }
  }
  
  #Checking that the date variable is in a date form
  if(!is.null(dateVar)){
    if(lubridate::is.Date(indData[, dateVar]) == FALSE &
      lubridate::is.POSIXt(indData[, dateVar]) == FALSE){
      stop(paste0(dateVar, " must be either a date or a date-time (POSIXt) object."))
    }
  }
  
  #Checking to make sure there is a date variable if ordered = TRUE
  if(ordered == TRUE & is.null(dateVar)){
    stop("If ordered = TRUE, then dateVar must be provided")
  }
  
  
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
    pairData2[, paste0(dateVar, ".Diff")] <- as.numeric(difftime(pairData2[, paste0(dateVar, ".2")],
                                              pairData2[, paste0(dateVar, ".1")],
                                              units = "days"))
    pairData2[, paste0(dateVar, ".DiffY")] <- pairData2[, paste0(dateVar, ".Diff")] / 365
  }
  
  if(ordered == TRUE){
    orderedData <- pairData2[pairData2[, paste0(dateVar, ".Diff")] > 0, ]
    return(orderedData)
  }else{
    return(pairData2)
  }
}


