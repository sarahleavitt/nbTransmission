
#' Calculates the effective reproductive number
#'
#' Uses probabilities to calculate the individual-level, time-level and overall average
#' effective reproductive numbers for an outbreak.
#' 
#' Add details section
#'
#' @param probs The name of the dateset with transmission probabilities
#' @param dateVar The variable name (in quotes) of the dates that the individuals are observed.
#' @param indIDVar The variable name (in quotes) of the individual ID.
#' @param pVar The variable name (in quotes) of the transmission probabilities.
#' @param timeFrame The time frame used to calculate Rt. One of "days", "weeks", "months", "years".
#' @param rangeForAvg A vector with the start and ending time period to be used to calculate the
#' average effective reproductive number.
#'
#' @return A list with three dataframes: one with the individual-level reproductive numbers,
#'  one with the time-level reproductive numbers, and one with the average effective reproductive
#'  number. ADD MORE DETAILS ON COLUMNS
#'
#' @examples
#' #Insert example here
#'
#' @import dplyr
#' 
#' @export

calcR <- function(probs, dateVar, indIDVar, pVar,
                   timeFrame = c("days", "months", "weeks", "years"),
                   rangeForAvg = NULL){
  
  #Creating correctly named variables
  probs <- as.data.frame(probs)
  probs$indID.1 <- probs[, paste0(indIDVar, ".1")]
  probs$indID.2 <- probs[, paste0(indIDVar, ".2")]
  probs$date.1 <- probs[, paste0(dateVar, ".1")]
  probs$date.2 <- probs[, paste0(dateVar, ".2")]
  probs$p <- probs[, pVar]
  
  #Calculating the individual-level reproductive number
  ri <- calcRi(probs)
  
  #Calculating the time-level reproductibe number
  rt <- calcRt(ri, timeFrame)
  
  #Calculating the average effective reproductive number
  rtAvg <- calcRtAvg(rt, rangeForAvg)
  

  return(list("RiDf" = ri, "RtDf" = rt[[1]], "RtAvg" = rtAvg))
}



#' Calculates Ri
#'
#' Uses probabilities to calculate the individual-level reproductive number
#'
#' @param probs The name of the dateset with transmission probabilities
#'
#' @return Individual-level dataframe with the following variables: indID, date, Ri,
#'      iden, nInfectors, nInfectees
#'
#' @import dplyr
#' 
#' @export


calcRi <- function(probs){
  
  #### Calculating identifiability of links ####
  
  #Calculating variance to get 1-q identifiability measure from Teunis et al.
  totalV <- (probs
             %>% mutate(var = p * (1 - p))
             %>% group_by(indID.2)
             %>% summarize(iden =  1- sum(var, na.rm = TRUE),
                           nInfectors = n(),
                           date = first(date.2))
             %>% rename(indID = indID.2)
  )
  
  #### Calculating reproductive number ####
  
  #Calculating the individual level reproductive number
  ri <- (probs
         %>% group_by(indID.1)
         %>% rename(indID = indID.1)
         %>% summarize(Ri = sum(p),
                       nInfectees = n(),
                       date = first(date.1))
         %>% full_join(totalV, by = c("indID", "date"))
         #If nInfectees is missing, this is the latest case
         #If nInfectors is missing, this is the earliest case
         #If Ri is missing, it is the last case so the value should be 0
         %>% tidyr::replace_na(list(nInfectees = 0, nInfectors = 0, Ri = 0))
         %>% select(indID, date, Ri, iden, nInfectors, nInfectees)
  )
  
  return(ri)
}



#' Calculates Rt
#'
#' Uses probabilities to calculate the time-level effective reproductive number from individual-level
#' reproductive numbers
#'
#' @param riData The name of the dateset with individual-level reproductive numbers.
#' @param timeFrame The time frame used to calculate Rt. One of "days", "weeks", "months", "years".
#'
#' @return Time-level dataframe with the following variables: timeRank = the rank of the time interval
#' used to calculate the reproductive number and Rt = the reproductive number for that time interval.
#'
#' @import dplyr
#' 
#' @export

calcRt <- function(riData, timeFrame = c("days", "weeks", "months", "years")){
  
  #### Creating ranks for the time interval ####
  
  #Finding the full range of time frames
  times <- seq(min(riData$date), max(riData$date), by = timeFrame)
  
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
    riData2 <- (riData 
               %>% mutate(time = format(date, timeFormat))
               %>% left_join(timeDf, by = "time")
    )
  }
  
  
  #I need to deal with weeks differently so that I can figure out the proper time frame
  if(timeFrame == "weeks"){
    #Creating a dataframe of week labels and their ranks
    timeDf <- cbind.data.frame(time = format(times, "%Y-%m-%d"),
                               timeRank = 1:length(times), stringsAsFactors = FALSE)
    
    #Ranking the days
    days <- seq(min(riData$date), max(riData$date), by = "days")
    dayDf <- cbind.data.frame(day = format(days, "%Y-%m-%d"),
                              dayRank = 1:length(days), stringsAsFactors = FALSE)
    
    #Finding the week rank based on the day rank / 7 and then merging back
    #with the time data frame to get the week label
    riData2 <- (riData 
                %>% mutate(day = format(date, "%Y-%m-%d"))
                %>% left_join(dayDf, by = "day")
                %>% mutate(timeRank = ceiling(dayRank / 7))
                %>% full_join(timeDf, by = "timeRank")
                %>% select(-day, -dayRank)
    )
  }
  
  
  #### Calculating Rt ####
  
  rt <- (riData2
         %>% group_by(timeRank)
         %>% summarize(Rt = mean(Ri))
         %>% ungroup()
  )
  
  return(list(rt))
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
    rtAvg <- mean(rtData$Rt)
    warning("Please choose the stable portion of the outbreak to calculate the average Rt")
  }else{
    rtCut <- rtData %>% filter(timeRank > rangeForAvg[1] & timeRank < rangeForAvg[2])
    rtAvg <- mean(rtCut$Rt)
  }
  
  return(as.data.frame(rtAvg))
}








