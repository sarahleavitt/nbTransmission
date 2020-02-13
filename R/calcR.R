
#' Estimates the effective reproductive number
#'
#' The function \code{estimateR} uses the relative transmission probabilities to estimate
#' the individual-level, time-level, and average effective reproductive numbers
#' for an outbreak.
#' 
#' The effective reproductive number is the average number of cases an infectious case
#' will produce in a population of both susceptible and non-susceptibe individuals.
#' The rational behind this reproductive number estimation is Wallinga and Teunis (2004)
#' where the individual-level reproductive number is estimated by summing
#' the relative probability that the individual infected any other individual. 
#' 
#' If \eqn{p_{ij}} equals the relative probability that case \eqn{i} infected
#' case \eqn{j}, then the individual-level reproductive number (\eqn{R_i}) is calculated by:
#'
#' \deqn{\sum_{m \ne i} {p_{im}}}
#' 
#' The time-level reproductive number is then estimated by averaging the individual-level
#' reproductive numbers for all individuals observed in the time frame (can specify days,
#' weeks, months, years).
#' 
#' Finally, the time-level reproductive numbers are averaged to
#' estimate the average effective reproductive number within \code{rangeForAvg}.
#' To get the best estimate of the average effective reproductive number, one should
#' only consider the stable portion of the outbreak (exclude the beginning and end).
#' 
#' If \code{bootSamples > 0}, bootstrap confidence intervals will be estimated for
#' both the time-level and average reproductive numbers using parametric bootstrapping.
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
#' \code{<dateVar>.2}) which must be date or date-time (POSIXt) objects.
#' @param pVar The column name (in quotes) of the transmission probabilities.
#' @param timeFrame The time frame used to calculate Rt
#' (one of \code{"days", "months", "weeks", "years"}).
#' @param rangeForAvg A vector with the start and ending time period to be used to calculate
#' the average effective reproductive number.
#' @param bootSamples The number of bootstrap samples; if 0, then no confidence intervals
#' are calculated.
#' @param alpha The alpha level for the confidence intervals.
#'
#' @return A list with five elements:
#'\enumerate{
#'   \item \code{RiDf} - a data frame with the individual-level reproductive numbers. Column names:
#'      \itemize{
#'        \item \code{<indIDVar>} - the individual ID with name specified.
#'        \item \code{<dateVar>} - the date the individual was observed with name specified.
#'        \item \code{Ri} - the individual-level reproductive number.
#'        \item \code{nInfectees} - the number of possible infectees for this individual.
#'      }
#'   \item \code{RtDf} - a data frame with the time-level reproductive numbers. Column names:
#'      \itemize{
#'        \item \code{time} - the time frame corresponding to the reproductive number estimate 
#'        (day for "days" and "weeks", month for "months", year for "years").
#'        \item \code{timeRank} - the rank of the time frame.
#'        \item \code{Rt} - the time-level reproductive number for this time frame.
#'        \item \code{ciLower} - lower bound of confidence interval for Rt
#'         (only if bootSamples > 0).
#'        \item \code{ciUpper} - upper bound of confidence interval for Rt 
#'        (only if bootSamples > 0).
#'      }
#'   \item \code{RtAvgDf} - a data frame with the average effective reproductive. Column names:
#'      \itemize{
#'        \item \code{RtAvg} - the average time-level reproductive number between the range
#'         specified in \code{rangeForAvg}.
#'        \item \code{ciLower} - lower bound of confidence interval for Rt
#'         (only if bootSamples > 0).
#'        \item \code{ciUpper} - upper bound of confidence interval for Rt 
#'        (only if bootSamples > 0).
#'      }
#'   \item \code{timeFrame} - a vector with the timeFrame input
#'   \item \code{rangeForAvg} - a vector with the rangeForAvg input
#' }
#' 
#' 
#' @examples
#' 
#' ## Use the nbResults data frame included in the package which has the results
#' ## of the nbProbabilities() function on a TB-like outbreak.
#' 
#' ## Getting initial estimates of the reproductive number
#' # (without specifying rangeForAvg and without confidence intervals)
#' rInitial <- estimateR(nbResults, dateVar = "infectionDate",
#'                indIDVar = "individualID", pVar = "pScaled",
#'                timeFrame = "months")
#'                
#' ## Finding the stable portion of the outbreak for rangeForAvg using plot
#' cut1 <- 25
#' cut2 <- 125
#' 
#' ## NOT RUN ##
#' # ggplot(data = rInitial$RtDf, aes(x = timeRank, y = Rt)) +
#' #   geom_point() +
#' #   geom_line() +
#' #   geom_hline(data = rInitial$RtAvgDf, aes(yintercept = RtAvg), size = 0.7) +
#' #   geom_vline(aes(xintercept = cut1), linetype = 2, size = 0.7) +
#' #   geom_vline(aes(xintercept = cut2), linetype = 2, size = 0.7)
#'   
#' ## Finding the final reproductive number estimates with confidence intervals
#' # NOTE should run with bootSamples > 10.
#' rFinal <- estimateR(nbResults, dateVar = "infectionDate",
#'              indIDVar = "individualID", pVar = "pScaled",
#'              timeFrame = "months", rangeForAvg = c(cut1, cut2),
#'              bootSamples = 2, alpha = 0.05)
#'
#' rFinal$RtAvgDf
#' 
#' @seealso \code{\link{nbProbabilities}} \code{\link{estimateRi}}
#'  \code{\link{estimateRt}} \code{\link{estimateRtAvg}}
#' 
#' 
#' @references 
#' Wallinga J, Teunis P. Different epidemic curves for severe acute respiratory
#' syndrome reveal similar impacts of control measures.
#' \emph{American Journal of Epidemiology}. 2004 Sep 15;160(6):509-16.
#' 
#' @export

estimateR <- function(df, indIDVar, dateVar, pVar,
                   timeFrame = c("days", "months", "weeks", "years"),
                   rangeForAvg = NULL, bootSamples = 0, alpha = 0.05){
  
  
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
  
  #Checking that the date variables are in a date form
  if(lubridate::is.Date(df[, dateVar1]) == FALSE &
     lubridate::is.POSIXt(df[, dateVar1]) == FALSE){
    stop(paste0(dateVar1, " must be either a date or a date-time (POSIXt) object."))
  }
  if(lubridate::is.Date(df[, dateVar2]) == FALSE &
     lubridate::is.POSIXt(df[, dateVar2]) == FALSE){
    stop(paste0(dateVar2, " must be either a date or a date-time (POSIXt) object."))
  }
  
  #Making sure timeFrame is one of the options
  if(!timeFrame %in% c("days", "months", "weeks", "years")){
    stop(paste0("timeFrame must be one of: ",
                paste0( c("days", "months", "weeks", "years"), collapse = ", ")))
  }

  
  
  #Calculating the individual-level reproductive number
  riEst <- estimateRi(df, pVar = pVar, indIDVar = indIDVar, dateVar = dateVar)
  
  #Calculating the time-level reproductibe number
  rtEst <- estimateRt(riEst, dateVar = dateVar, timeFrame = timeFrame)
  
  #Calculating the average effective reproductive number
  rtAvgEst <- estimateRtAvg(rtEst, rangeForAvg)
  
  if(bootSamples == 0){
    return(list("RiDf" = riEst, "RtDf" = rtEst, "RtAvgDf" = rtAvgEst,
                "timeFrame" = timeFrame, "rangeForAvg" = rangeForAvg))}
  
  if (bootSamples > 0){
    
    ## Finding the CI for Rt ##
    
    #Use manual list instead if replicate for progress bar
    pb <- utils::txtProgressBar(min = 0, max = bootSamples, style = 3)
    bootRt <- data.frame(c("timeRank" = integer(), "Rt" = numeric(),
                           "time" = character()), "rep" = integer(),
                         stringsAsFactors = FALSE)
    
    for(i in 1:bootSamples){
      oneRep <- estimateRt(simulateRi(df, riEst, pVar = pVar, indIDVar = indIDVar),
                       dateVar = dateVar, timeFrame = timeFrame)
      oneRep$rep <- i
      bootRt <- rbind(bootRt, oneRep)
      utils::setTxtProgressBar(pb, i)
    }
    
    #Finding the quantiles of the distribution of Rt
    cilb <- stats::aggregate(bootRt$Rt, by = list(bootRt$timeRank),
                      stats::quantile, probs = 1-alpha/2)
    names(cilb) <- c("timeRank", "lb")
    ciub <- stats::aggregate(bootRt$Rt, by = list(bootRt$timeRank),
                      stats::quantile, probs = alpha/2)
    names(ciub) <- c("timeRank", "ub")
    
    #Combing the quantiles with Rt data
    ciRt1 <- merge(rtEst, cilb, by = "timeRank", all = TRUE)
    ciRt2 <- merge(ciRt1, ciub, by = "timeRank", all = TRUE)
    #Calculating the CI bounds
    ciRt2$ciLower <- ifelse(ciRt2$Rt - (ciRt2$lb - ciRt2$Rt) > 0,
                            ciRt2$Rt - (ciRt2$lb - ciRt2$Rt), 0)
    ciRt2$ciUpper <- ciRt2$Rt - (ciRt2$ub - ciRt2$Rt)
    ciDataRt <- ciRt2[, !names(ciRt2) %in% c("lb", "ub")]
    
    #Reordering columns
    ciDataRt <- ciDataRt[, c("time", "timeRank", "Rt", "ciLower", "ciUpper")]
    
    
    ## Finding the CI for RtAvg ##
    
    #Calculating RtAvg for each bootstrap Rt sample
    bootRtAvgL <- by(bootRt, INDICES = list(bootRt$rep),
                    FUN = estimateRtAvg, rangeForAvg = rangeForAvg)
    bootRtAvg <- do.call(c, bootRtAvgL)
    
    #Calculating the CI values
    ciLower <- rtAvgEst - (stats::quantile(bootRtAvg, 1-alpha/2) - rtAvgEst)
    ciUpper <- rtAvgEst - (stats::quantile(bootRtAvg, alpha/2) - rtAvgEst)
    ciDataRtAvg <- cbind.data.frame(rtAvgEst, ciLower, ciUpper, row.names = NULL)
    names(ciDataRtAvg) <- c("RtAvg", "ciLower", "ciUpper")
    
    return(list("RiDf" = riEst, "RtDf" = ciDataRt, "RtAvgDf" = ciDataRtAvg,
                "timeFrame" = timeFrame, "rangeForAvg" = rangeForAvg))
  }
}



#' Estimates individual-level reproductive numbers
#'
#' The function \code{estimateRi} uses relative transmission probabilities to estimate the
#' individual-level reproductive number.
#' 
#' This function is meant to be called by \code{\link{estimateR}}
#' which estimates the individual-level, time-level, and average reproductive numbers, 
#' but it can also be run directly.
#'
#' @param df The name of the dateset with transmission probabilities
#' @param indIDVar The variable name (in quotes) of the individual ID varaibles 
#' (data frame \code{df} must have variables called \code{<indIDVar>.1}
#' and \code{<indIDVar>.2}).
#' @param dateVar The variable name (in quotes) of the dates that the individuals
#' are observed (data frame \code{df}  must have variables called \code{<dateVar>.1}
#' and \code{<dateVar>.2}).
#' @param pVar The variable name (in quotes) of the transmission probabilities.
#'
#' @return A data frame with the individual-level reproductive numbers. Column names:
#'      \itemize{
#'        \item \code{<indIDVar>} - the individual ID with name specified.
#'        \item \code{<dateVar>} - the date the individual was observed with name specified.
#'        \item \code{Ri} - the individual-level reproductive number.
#'        \item \code{nInfectees} - the number of possible infectees for this individual.
#'      }
#'      
#' @seealso \code{\link{estimateR}} \code{\link{estimateRt}} \code{\link{estimateRtAvg}}
#' 
#' @export


estimateRi <- function(df, indIDVar, dateVar, pVar){
  
  df <- as.data.frame(df)
  #Creating variables with the individual indID and date variables
  indIDVar1 <- paste0(indIDVar, ".1")
  indIDVar2 <- paste0(indIDVar, ".2")
  dateVar1 <- paste0(dateVar, ".1")
  dateVar2 <-  paste0(dateVar, ".2")
  
  #### Calculating individual-level reproductive number ####
  
  #Summing the scaled probabilties to get individual-level reproductive numer
  sumP <- stats::aggregate(df[, pVar], by = list(df[, indIDVar1]), sum)
  names(sumP) <- c(indIDVar, "Ri")
  #Finding the number of possible infectors per infectee
  nInf <- as.data.frame(table(df[, indIDVar1]))
  names(nInf) <- c(indIDVar, "nInfectees")
  #Extracting the date of infection for each infectee
  dateInfo <- df[!duplicated(df[, indIDVar1]), c(indIDVar1, dateVar1)]
  names(dateInfo) <- c(indIDVar, dateVar)
  
  #Combining these variables to create Ri data frame
  ri1 <- merge(dateInfo, sumP, by = indIDVar, all = TRUE)
  ri <- merge(ri1, nInf, by = indIDVar, all = TRUE)
  #If nInfectees is missing and RI is missing, this is the latest case
  #so the values should be 0
  ri[is.na(ri$Ri), "Ri"] <- 0
  ri[is.na(ri$nInfectees), "nInfectees"] <- 0
  
  return(ri)
}



#' Estimates time-level reproductive numbers
#'
#' The function \code{estimateRt} estimates the time-level effective reproductive number
#' from individual-level reproductive numbers.
#' 
#' This function is meant to be called by \code{\link{estimateR}}
#' which estimates the individual-level and time-level, and average reproductive numbers, 
#' but it can also be run directly.
#'
#' @param riData The name of the dateset with individual-level reproductive numbers.
#' @param dateVar The variable name (in quotes) of the dates that the individuals are observed
#' (data frame \code{riData}  must have a variable called \code{<dateVar>}).
#' @param timeFrame The time frame used to calculate Rt.
#'
#' @return A data frame with the time-level reproductive numbers. Column names:
#'      \itemize{
#'        \item \code{time} - the time frame corresponding to the reproductive number estimate 
#'        (day for "days" and "weeks", month for "months", year for "years").
#'        \item \code{timeRank} - the rank of the time frame.
#'        \item \code{Rt} - the time-level reproductive number for this time frame.
#'      }
#' 
#' @seealso \code{\link{estimateR}} \code{\link{estimateRi}} \code{\link{estimateRtAvg}}
#' 
#' @export

estimateRt <- function(riData, dateVar, timeFrame = c("days", "weeks", "months", "years")){
  
  #### Creating ranks for the time interval ####
  
  #Finding the full range of time frames
  times <- seq(min(riData[, dateVar]), max(riData[, dateVar]), by = timeFrame)
  
  #Defining the format for the time frame
  if(timeFrame == "days"){
    timeFormat <- "%Y-%m-%d"
  }else if(timeFrame == "months"){
    timeFormat <- "%Y-%m"
  }else if(timeFrame == "years"){
    timeFormat <- "%Y"
  }
  
  if(timeFrame %in% c("days", "months", "years")){
    timeDf <- cbind.data.frame(time = format(times, timeFormat),
                               timeRank = 1:length(times), stringsAsFactors = FALSE)
    riData$time <- format(riData[, dateVar], timeFormat)
    riDataFinal <- merge(riData, timeDf, by = "time", all.x = TRUE)
  }
  
  #Need to deal with weeks differently to figure out the proper time frame
  if(timeFrame == "weeks"){
    #Creating a data frame of week labels and their ranks
    timeDf <- cbind.data.frame(time = format(times, "%Y-%m-%d"),
                               timeRank = 1:length(times), stringsAsFactors = FALSE)
    
    #Ranking the days
    days <- seq(min(riData[, dateVar]), max(riData[, dateVar]), by = "days")
    dayDf <- cbind.data.frame(day = format(days, "%Y-%m-%d"),
                              dayRank = 1:length(days), stringsAsFactors = FALSE)
    
    #Finding the week rank based on the day rank / 7 and then merging back
    #with the time data frame to get the week label
    riData$day <- format(riData[, dateVar], "%Y-%m-%d")
    riData2 <- merge(riData, dayDf, by = "day", all.x = TRUE)
    riData2$timeRank <- ceiling(riData2$dayRank / 7)
    riData3 <- merge(riData2, timeDf, by = "timeRank")
    riDataFinal <- riData3[, !names(riData3) %in% c("day", "dayRank")]
  }
  
  ## Calculating Rt ##
  
  #Finding the mean Ri value by time period (Rt)
  meanRi <- stats::aggregate(riDataFinal$Ri, by = list(riDataFinal$timeRank), mean)
  names(meanRi) <- c("timeRank", "Rt")
  #Extracting the time for each Rt value
  timeInfo <- riDataFinal[!duplicated(riDataFinal$timeRank), c("timeRank", "time")]
  #Combining these data frames to create Rt data frame
  rt <- merge(meanRi, timeInfo, by = "timeRank", all = TRUE)
  
  #Ordering columns
  rt <- rt[, c("time", "timeRank", "Rt")]
  
  return(rt)
}



#' Estimates the average effective reproductive number
#'
#' Averages the time-level reproductive numbers within a certain range to estimate the overall
#' reproductive number for an oubreak.
#' 
#' This function is meant to be called by \code{\link{estimateR}}
#' which estimates the individual-level and time-level, and average reproductive numbers,
#' but it can also be run directly.
#'
#' @param rtData The name of the dateset with time-level reproductive numbers.
#' @param rangeForAvg A vector with the start and ending time period to be used to calculate the
#' average effective reproductive number.
#'
#' @return A data frame with the average effective reproductive. Column names:
#'      \itemize{
#'        \item \code{RtAvg} - the average time-level reproductive number between the range
#'         specified in \code{rangeForAvg}.
#'      }
#' 
#' @seealso \code{\link{estimateR}} \code{\link{estimateRi}} \code{\link{estimateRt}}
#'  
#' @export


estimateRtAvg <- function(rtData, rangeForAvg = NULL){
  
  #Printing a note if you do not specify a range for the average 
  if(is.null(rangeForAvg)){
    #Finding the mean of the time-level reproductive numbers
    RtAvg <- mean(rtData$Rt)
    message("Please choose the stable portion of the outbreak to calculate the average Rt")
  }else{
    #Restricting to the time-level reproductive numbers within the rangeForAverage.
    rtCut <- rtData[rtData$timeRank > rangeForAvg[1] &
                      rtData$timeRank < rangeForAvg[2], ]
    RtAvg <- mean(rtCut$Rt, na.rm = TRUE)
  }
  
  return(as.data.frame(RtAvg))
}


#Function that simulates Ri from the probabilities
simulateRi <- function(df, riEst, pVar, indIDVar){
  
  #Creating variables with the individual indID and date variables
  indIDVar1 <- paste0(indIDVar, ".1")
  indIDVar2 <- paste0(indIDVar, ".2")
  
  #Making a matrix of probabilities for speedy calculation
  dfM <- unclass(stats::xtabs(df[, pVar] ~ df[, indIDVar1] + df[, indIDVar2]))
  Ri <- apply(dfM, 1, function(x) poisbinom::rpoisbinom(n = 1, pp = x))
  riNew <- cbind.data.frame(names(Ri), Ri, stringsAsFactors = "FALSE")
  names(riNew) <- c(indIDVar, "Ri")
  
  #Merging the new Ri value with the individual information
  riEst$indIDChar <- as.character(riEst[, indIDVar])
  riEst2 <- riEst[!names(riEst) %in% c("Ri")]
  riNew2 <- merge(riEst2, riNew, by = c("indIDChar" = indIDVar))
  riNew2[is.na(riNew2$Ri), "Ri"] <- 0
  riNew3 <- riNew2[, !names(riNew2) %in% c("indIDChar")]
  
  return(riNew3)
}




