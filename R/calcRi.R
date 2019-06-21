
#' Calculates Ri
#'
#' Uses probabilities to calculate the individual level reproductive number
#'
#' @param probs dateset with transmission probabilities
#' @param dateVar the variable name (in quotes) of the date that the individual is observed
#' @param indIDVar the variable name (in quotes) of the id variable
#' @param pVar the variable name (in quotes) of the transmission probabilities
#'
#' @return Individual-level dataframe with the following variables: indID, Ri, date,
#'      iden, nInfectors, month, year, dayR, monthR, yearR
#'
#' @examples
#' #Insert example here
#'
#' @export


calcRi <- function(probs, dateVar, indIDVar, pVar){
  
  #Creating correctly named variables
  probs <- as.data.frame(probs)
  probs$date.1 <- probs[, paste0(dateVar, ".1")]
  probs$date.2 <- probs[, paste0(dateVar, ".2")]
  probs$timeDiff <- as.numeric(difftime(probs$date.2, probs$date.1, units = "days"))
  probs$indID.1 <- probs[, paste0(indIDVar, ".1")]
  probs$indID.2 <- probs[, paste0(indIDVar, ".2")]
  probs$p <- probs[, pVar]
  
  
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
         %>% mutate(label = label)
         #If nInfectees is missing, this is the latest case
         #If nInfectors is missing, this is the earliest case
         #If Ri is missing, it is the last case so the value should be 0
         %>% replace_na(list(nInfectees = 0, nInfectors = 0, Ri = 0))
  )
  
  
  #### Creating day, month, and year variables ####
  
  #Finding all days, months, and years represented in the outbreak
  days <- seq(min(ri$date), max(ri$date), by = "days")
  daysDf <- cbind.data.frame(date = format(days, "%Y-%m-%d"),
                             dayR = 1:length(days), stringsAsFactors = FALSE)
  
  months <- format(seq(min(ri$date), max(ri$date), by = "months"), "%Y-%m")
  monthsDf <- cbind.data.frame(month = months, monthR = 1:length(months), stringsAsFactors = FALSE)
  
  years <- year(seq(min(ri$date), max(ri$date), by = "years"))
  yearsDf <- cbind.data.frame(year = years, yearR = 1:length(years), stringsAsFactors = FALSE)

  
  #Combining the year and month rankings with the raw data
  riFull <- (ri
             %>% mutate(month = format(date, "%Y-%m"),
                        year = year(date),
                        date = format(date, "%Y-%m-%d"))
             %>% left_join(daysDf, by = "date")
             %>% left_join(monthsDf, by = "month")
             %>% left_join(yearsDf, by = "year")
  )
  
  return(riFull)
}


