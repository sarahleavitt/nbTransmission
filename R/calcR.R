
#' Calculates the effective reproductive number
#'
#' The function \code{calcR} uses the relative transmission probabilities to estimate
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
#' both the time-level and average reproductive numbers.
#' 
#' 
#' @param probs The name of the dateset with transmission probabilities
#' @param indIDVar The variable name (in quotes) of the individual ID varaibles
#' (dataframe \code{probs} 
#' must have variables called \code{<indIDVar>.1} and \code{<indIDVar>.2}).
#' @param dateVar The variable name (in quotes) of the dates that the individuals are
#' observed (dataframe \code{probs}  must have variables called \code{<dateVar>.1} and
#' \code{<dateVar>.2}).
#' @param pVar The variable name (in quotes) of the transmission probabilities.
#' @param timeFrame The time frame used to calculate Rt.
#' @param rangeForAvg A vector with the start and ending time period to be used to calculate
#' the average effective reproductive number.
#' @param bootSamples The number of bootstrap samples; if 0, then no confidence intervals
#' are calculated.
#' @param alpha The alpha level for the confidence intervals.
#'
#' @return A list with three dataframes:
#'\enumerate{
#'   \item \code{RiDf} - a dataframe with the individual-level reproductive numbers. Column names:
#'      \itemize{
#'        \item \code{<indIDVar>} - the individual ID with name specified.
#'        \item \code{<dateVar>} - the date the individual was observed with name specified.
#'        \item \code{Ri} - the individual-level reproductive number.
#'        \item \code{nInfectees} - the number of possible infectees for this individual.
#'      }
#'   \item \code{RtDf} - a dataframe with the time-level reproductive numbers. Column names:
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
#'   \item \code{RtAvgDf} - a dataframe with the average effective reproductive. Column names:
#'      \itemize{
#'        \item \code{RtAvg} - the average time-level reproductive number between the range
#'         specified in \code{rangeForAvg}.
#'        \item \code{ciLower} - lower bound of confidence interval for Rt
#'         (only if bootSamples > 0).
#'        \item \code{ciUpper} - upper bound of confidence interval for Rt 
#'        (only if bootSamples > 0).
#'      }
#' }
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
#' #NOTE should run with nReps > 1
#' covariates = c("Z1", "Z2", "Z3", "Z4", "timeCat")
#' resGen <- calcProbabilities(orderedPair = orderedPair,
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
#' ## Getting initial estimates of the reproductive number
#' # (ithout specifying rangeForAvg and without confidence intervals)
#' rInitial <- calcR(allProbs, dateVar = "infectionDate",
#'                indIDVar = "individualID", pVar = "pScaled",
#'                timeFrame = "months")
#'                
#' ## Finding the stable portion of the outbreak for rangeForAvg
#' rt <- rInitial$RtDf
#' totalTime <- max(rt$timeRank, na.rm = TRUE) - min(rt$timeRank, na.rm = TRUE)
#' monthCut1 <- ceiling(0.15 * totalTime)
#' monthCut2 <- ceiling(0.85 * totalTime)
#' 
#' ## NOT RUN ##
#' # ggplot(data = rt, aes(x = timeRank, y = Rt)) +
#' #   geom_point() +
#' #   geom_line() +
#' #   geom_hline(data = rInitial$RtAvgDf, aes(yintercept = RtAvg), size = 0.7) +
#' #   geom_vline(aes(xintercept = monthCut1), linetype = 2, size = 0.7) +
#' #   geom_vline(aes(xintercept = monthCut2), linetype = 2, size = 0.7)
#'   
#' ## Finding the final reproductive number estimates with confidence intervals
#' # NOTE should run with bootSamples > 10.
#' rFinal <- calcR(allProbs, dateVar = "infectionDate",
#'              indIDVar = "individualID", pVar = "pScaled",
#'              timeFrame = "months",
#'              rangeForAvg = c(monthCut1, monthCut2),
#'              bootSamples = 10, alpha = 0.05)

#' rFinal$RtAvgDf
#'   
#'   
#' 
#' @references 
#' Wallinga J, Teunis P. Different epidemic curves for severe acute respiratory syndrome reveal
#' similar impacts of control measures. American Journal of epidemiology.
#' 2004 Sep 15;160(6):509-16.
#' 
#' @export

calcR <- function(probs, indIDVar, dateVar, pVar,
                   timeFrame = c("days", "months", "weeks", "years"),
                   rangeForAvg = NULL, bootSamples = 0, alpha = 0.05){
  
  probs <- as.data.frame(probs)
  #Creating variables with the individual indID and date variables
  indIDVar1 <- paste0(indIDVar, ".1")
  indIDVar2 <- paste0(indIDVar, ".2")
  dateVar1 <- paste0(dateVar, ".1")
  dateVar2 <-  paste0(dateVar, ".2")
  
  #Calculating the individual-level reproductive number
  riEst <- calcRi(probs, pVar = pVar, indIDVar = indIDVar, dateVar = dateVar)
  
  #Calculating the time-level reproductibe number
  rtEst <- calcRt(riEst, dateVar = dateVar, timeFrame = timeFrame)
  
  #Calculating the average effective reproductive number
  rtAvgEst <- calcRtAvg(rtEst, rangeForAvg)
  
  if(bootSamples == 0){
    return(list("RiDf" = riEst, "RtDf" = rtEst, "RtAvgDf" = rtAvgEst))}
  
  if (bootSamples > 0){
    
    ## Finding the CI for Rt ##
    
    #Use manual list instead if replicate for progress bar
    pb <- utils::txtProgressBar(min = 0, max = bootSamples, style = 3)
    bootRt <- data.frame(c("timeRank" = integer(), "Rt" = numeric(),
                           "time" = character()), "rep" = integer())
    for(i in 1:bootSamples){
      oneRep <- calcRt(simulateRi(probs, riEst, pVar = pVar, indIDVar = indIDVar),
                       dateVar = dateVar, timeFrame = "months")
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
                    FUN = calcRtAvg, rangeForAvg = rangeForAvg)
    bootRtAvg <- do.call(c, bootRtAvgL)
    
    #Calculating the CI values
    ciLower <- rtAvgEst - (stats::quantile(bootRtAvg, 1-alpha/2) - rtAvgEst)
    ciUpper <- rtAvgEst - (stats::quantile(bootRtAvg, alpha/2) - rtAvgEst)
    ciDataRtAvg <- cbind.data.frame(rtAvgEst, ciLower, ciUpper, row.names = NULL)
    names(ciDataRtAvg) <- c("RtAvg", "ciLower", "ciUpper")
    
    return(list("RiDf" = riEst, "RtDf" = ciDataRt, "RtAvgDf" = ciDataRtAvg))
  }
}



#' Estimates individual-level reproductive numbers
#'
#' The function \code{calcRi} uses relative transmission probabilities to estimate the
#' individual-level reproductive number.
#' 
#' This function is meant to be called by \code{\link{calcR}}
#' which estimates the individual-level, time-level, and average reproductive numbers, 
#' but it can also be run directly.
#'
#' @param probs The name of the dateset with transmission probabilities
#' @param indIDVar The variable name (in quotes) of the individual ID varaibles 
#' (dataframe \code{probs} must have variables called \code{<indIDVar>.1}
#' and \code{<indIDVar>.2}).
#' @param dateVar The variable name (in quotes) of the dates that the individuals
#' are observed (dataframe \code{probs}  must have variables called \code{<dateVar>.1}
#' and \code{<dateVar>.2}).
#' @param pVar The variable name (in quotes) of the transmission probabilities.
#'
#' @return A dataframe with the individual-level reproductive numbers. Column names:
#'      \itemize{
#'        \item \code{<indIDVar>} - the individual ID with name specified.
#'        \item \code{<dateVar>} - the date the individual was observed with name specified.
#'        \item \code{Ri} - the individual-level reproductive number.
#'        \item \code{nInfectees} - the number of possible infectees for this individual.
#'      }
#'      
#' @seealso \code{\link{calcR}}
#' 
#' @export


calcRi <- function(probs, indIDVar, dateVar, pVar){
  
  probs <- as.data.frame(probs)
  #Creating variables with the individual indID and date variables
  indIDVar1 <- paste0(indIDVar, ".1")
  indIDVar2 <- paste0(indIDVar, ".2")
  dateVar1 <- paste0(dateVar, ".1")
  dateVar2 <-  paste0(dateVar, ".2")
  
  #### Calculating individual-level reproductive number ####
  
  #Summing the scaled probabilties to get individual-level reproductive numer
  sumP <- stats::aggregate(probs[, pVar], by = list(probs[, indIDVar1]), sum)
  names(sumP) <- c(indIDVar, "Ri")
  #Finding the number of possible infectors per infectee
  nInf <- as.data.frame(table(probs[, indIDVar1]))
  names(nInf) <- c(indIDVar, "nInfectees")
  #Extracting the date of infection for each infectee
  dateInfo <- probs[!duplicated(probs[, indIDVar1]), c(indIDVar1, dateVar1)]
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
#' The function \code{calcRt} estimates the time-level effective reproductive number
#' from individual-level reproductive numbers.
#' 
#' This function is meant to be called by \code{\link{calcR}}
#' which estimates the individual-level and time-level, and average reproductive numbers, 
#' but it can also be run directly.
#'
#' @param riData The name of the dateset with individual-level reproductive numbers.
#' @param dateVar The variable name (in quotes) of the dates that the individuals are observed
#' (dataframe \code{riData}  must have a variable called \code{<dateVar>}).
#' @param timeFrame The time frame used to calculate Rt.
#'
#' @return A dataframe with the time-level reproductive numbers. Column names:
#'      \itemize{
#'        \item \code{time} - the time frame corresponding to the reproductive number estimate 
#'        (day for "days" and "weeks", month for "months", year for "years").
#'        \item \code{timeRank} - the rank of the time frame.
#'        \item \code{Rt} - the time-level reproductive number for this time frame.
#'      }
#' 
#' @seealso \code{\link{calcR}}
#' 
#' @export

calcRt <- function(riData, dateVar, timeFrame = c("days", "weeks", "months", "years")){
  
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
  
  #I need to deal with weeks differently so that I can figure out the proper time frame
  if(timeFrame == "weeks"){
    #Creating a dataframe of week labels and their ranks
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
  #Combining these dataframes to create Rt data frame
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
#' This function is meant to be called by \code{\link{calcR}}
#' which estimates the individual-level and time-level, and average reproductive numbers,
#' but it can also be run directly.
#'
#' @param rtData The name of the dateset with time-level reproductive numbers.
#' @param rangeForAvg A vector with the start and ending time period to be used to calculate the
#' average effective reproductive number.
#'
#' @return A dataframe with the average effective reproductive. Column names:
#'      \itemize{
#'        \item \code{RtAvg} - the average time-level reproductive number between the range
#'         specified in \code{rangeForAvg}.
#'      }
#' 
#' @seealso \code{\link{calcR}}
#'  
#' @export


calcRtAvg <- function(rtData, rangeForAvg = NULL){
  
  #Printing a note if you do not specify a range for the average 
  if(is.null(rangeForAvg)){
    #Finding the mean of the time-level reproductive numbers
    RtAvg <- mean(rtData$Rt)
    print("Please choose the stable portion of the outbreak to calculate the average Rt")
  }else{
    #Restricting to the time-level reproductive numbers within the rangeForAverage.
    rtCut <- rtData[rtData$timeRank > rangeForAvg[1] &
                      rtData$timeRank < rangeForAvg[2], ]
    RtAvg <- mean(rtCut$Rt)
  }
  
  return(as.data.frame(RtAvg))
}


#Function that simulates Ri from the probabilities
simulateRi <- function(probs, riEst, pVar, indIDVar){
  
  #Creating variables with the individual indID and date variables
  indIDVar1 <- paste0(indIDVar, ".1")
  indIDVar2 <- paste0(indIDVar, ".2")
  
  #Making a matrix of probabilities for speedy calculation
  probsM <- unclass(stats::xtabs(probs[, pVar] ~ probs[, indIDVar1] + probs[, indIDVar2]))
  Ri <- apply(probsM, 1, function(x) poisbinom::rpoisbinom(n = 1, pp = x))
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




