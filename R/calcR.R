
#' Calculates the effective reproductive number
#'
#' Uses relative transmission probabilities to calculate the individual-level,
#' time-level and overall average effective reproductive numbers for an outbreak.
#' 
#' Add details section
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
#' @return A list with three dataframes: one with the individual-level reproductive numbers,
#'  one with the time-level reproductive numbers, and one with the average effective reproductive
#'  number. ADD MORE DETAILS ON COLUMNS
#'
#' @examples
#' #Insert example here
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



#' Calculates Ri
#'
#' Uses probabilities to calculate the individual-level reproductive number
#'
#' @param probs The name of the dateset with transmission probabilities
#' @param indIDVar The variable name (in quotes) of the individual ID varaibles 
#' (dataframe \code{probs} must have variables called \code{<indIDVar>.1} and \code{<indIDVar>.2}).
#' @param dateVar The variable name (in quotes) of the dates that the individuals are observed
#' (dataframe \code{probs}  must have variables called \code{<dateVar>.1} and \code{<dateVar>.2}).
#' @param pVar The variable name (in quotes) of the transmission probabilities.
#'
#' @return Individual-level dataframe with the following variables: indID, date, Ri,
#'      iden, nInfectors, nInfectees
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



#' Calculates Rt
#'
#' Uses probabilities to calculate the time-level effective reproductive number from individual-level
#' reproductive numbers
#'
#' @param riData The name of the dateset with individual-level reproductive numbers.
#' @param dateVar The variable name (in quotes) of the dates that the individuals are observed
#' (dataframe \code{riData}  must have a variable called \code{<dateVar>}).
#' @param timeFrame The time frame used to calculate Rt.
#'
#' @return Time-level dataframe with the following variables: timeRank = the rank of the time interval
#' used to calculate the reproductive number and Rt = the reproductive number for that time interval.
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
  
  return(rt)
}



#' Calculates Average Rt
#'
#' Uses probabilities to calculate the average time-level effective reproductive number.
#'
#' @param rtData The name of the dateset with time-level reproductive numbers.
#' @param rangeForAvg A vector with the start and ending time period to be used to calculate the
#' average effective reproductive number.
#'
#' @return The value of the average time-level reproductive number for the outbreak within
#' "rangeForAvg"
#'
#' @import dplyr
#' 
#' @export


calcRtAvg <- function(rtData, rangeForAvg = NULL){
  
  if(is.null(rangeForAvg)){
    RtAvg <- mean(rtData$Rt)
    warning("Please choose the stable portion of the outbreak to calculate the average Rt")
  }else{
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




