
#' Calculates the bootstrap confidence intervals for the effective reproductive number
#'
#' Uses probabilities to resample the individual-level reproductive numbers in order to
#' estimate bootstrap confidence intervals for the time-level and average effective
#' reproductive numbers.
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

bootstrapR <- function(probs, dateVar, indIDVar, pVar,
                       timeFrame = c("days", "months", "weeks", "years"),
                       rangeForAvg = NULL, B = 1000, alpha = 0.05){
  
  #Creating correctly named variables
  probs <- as.data.frame(probs)
  probs$indID.1 <- probs[, paste0(indIDVar, ".1")]
  probs$indID.2 <- probs[, paste0(indIDVar, ".2")]
  probs$date.1 <- probs[, paste0(dateVar, ".1")]
  probs$date.2 <- probs[, paste0(dateVar, ".2")]
  probs$p <- probs[, pVar]
  
  #Finding the point estimates
  riEst <- calcRi(probs)
  rtEst <- calcRt(riEst, timeFrame)[[1]]
  rtAvgEst <- calcRtAvg(rtEst, rangeForAvg)$rtAvg
  
  
  ## Finding the CI for Rt ##
  
  bootListRt <- replicate(B, calcRt(simulateRi(probs, riEst), timeFrame = "months"))
  bootRt <- bind_rows(bootListRt)
  
  ciDataRt <- (bootRt
               %>% group_by(timeRank)
               %>% summarize(lb = quantile(Rt, 1-alpha/2),
                             ub = quantile(Rt, alpha/2))
               %>% full_join(rtEst, by = "timeRank")
               %>% mutate(ciLower = ifelse(Rt - (lb - Rt) > 0, 
                                         Rt - (lb - Rt), 0),
                        ciUpper = Rt - (ub - Rt))
               %>% select(-lb, -ub)
  )
  
  
  ## Finding the CI for RtAvg ##
  
  #Calculating RtAvg for each bootstrap Rt sample
  bootRtAvg <- (bootRt
                %>% group_by(timeRank)
                %>% mutate(rep = 1:n())
                %>% group_by(rep)
                %>% do(calcRtAvg(., rangeForAvg))
                %>% pull(rtAvg)
  )
  
  ciLower <- rtAvgEst - (quantile(bootRtAvg, 1-alpha/2) - rtAvgEst)
  ciUpper <- rtAvgEst - (quantile(bootRtAvg, alpha/2) - rtAvgEst)
  ciDataRtAvg <- cbind.data.frame(RtAvg = rtAvgEst, ciLower = ciLower, 
                                  ciUpper = ciUpper, row.names = NULL)
  
  return(list("RiDf" = riEst, "RtCI" = ciDataRt, "RtAvgCI" = ciDataRtAvg))
}


#Function that simulates Ri from the probabilities
simulateRi <- function(probs, riEst){

  #Making a matrix of probabilities for speedy calculation
  probsM <- unclass(xtabs(probs$p ~ probs$indID.1 + probs$indID.2))
  Ri <- apply(probsM, 1, function(x) poisbinom::rpoisbinom(n = 1, pp = x))
  riNew <- cbind.data.frame(indID = names(Ri), Ri, stringsAsFactors = "FALSE")
  
  #Merging this data back with the original Ri values
  riNew2 <- (riEst
             %>% mutate(indIDChar = as.character(indID))
             %>% select(-Ri)
             %>% full_join(riNew, by = c("indIDChar" = "indID"))
             %>% replace_na(list(Ri = 0))
             %>% select(-indIDChar)
  )
  
  return(riNew2)
}


