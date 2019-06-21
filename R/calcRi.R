#Sarah Van Ness
#Boston University
#Dissertation

################################################################################
# This program creates a function that calculates the individual
# reproductive numbers
################################################################################


################## Function that calculated Ri ####################


#### Inputs ####

#label = label for run (represent scheme, method, goldStandard, and ratio)
#results = dataset with probabilities for each pair
#pVar = variable that represents the probabilities to assess (default is "pScaled")
#label = label of the run you want to analyze (default is NA)
#labelVar = variable that has the runIDs (default is "label")


#### Outputs ####

#Individual-level dataframe with the following variables:
  #individualID
  #Ri - the individual-level reproductive number
  #nInfectees - the number of possible infectees for this case (number of cases observed earlier)
  #outbreakID - the outbreakID for this case
  #observationDate - the date of observation of the case
  #iden - the identification value (how sure we are of this cases's infector) from Teunis et al. 
  #nInfectors - the number of possible infectors for this case (number of cases observed later)
  #label - the run label for this analysis (or missing if not applicable)
  #month - the year and month of observation
  #year - the year of observation
  #dayR - the rank of the day of observation
  #monthR - the rank of the year and month of observation
  #yearR - the rank of the year of observation



calcRepNum <- function(label = NA, results, labelVar = "label", pVar = "pScaled"){
  
  results <- as.data.frame(results)
  results$p <- results[, pVar]
  
  #If there is a label, subsetting to just that run
  if(!is.na(label) & !is.na(labelVar)){
    allM <- results[results[, labelVar] == label, ]
  }else{
    allM <- results
  }

  
  #### Calculating identifiability of links ####
  
  #If outbreakID is missing, setting it to 1 (only need this variable for simulations)
  if(!"outbreakID.2" %in% names(allM)){
    allM$outbreakID.1 = 1
    allM$outbreakID.2 = 1
  }
  #If runID is missing, setting it to 1_1 (only need this variable for simulations)
  if(!"runID" %in% names(allM)){
    allM$runID = "1_1"
  }
  
  #Calculating variance to get 1-q identifiability measure from Teunis et al.
  totalV <- (allM
             %>% mutate(var = p * (1 - p))
             %>% group_by(individualID.2)
             %>% summarize(iden =  1- sum(var, na.rm = TRUE),
                           nInfectors = n(),
                           outbreakID = first(outbreakID.2),
                           observationDate = first(observationDate.2))
             %>% rename(individualID = individualID.2)
  )
  
  
  #### Calculating reproductive number ####
  
  #Calculating the individual level reproductive number
  repNumI <- (allM
              %>% group_by(individualID.1)
              %>% rename(individualID = individualID.1)
              %>% summarize(Ri = sum(p),
                            nInfectees = n(),
                            outbreakID = first(outbreakID.1),
                            observationDate = first(observationDate.1),
                            runID = first(runID))
              %>% full_join(totalV, by = c("individualID", "outbreakID", "observationDate"))
              %>% mutate(label = label)
              #If nInfectees is missing, this is the latest case
              #If nInfectors is missing, this is the earliest case
              #If Ri is missing, it is the last case so the value should be 0
              %>% replace_na(list(nInfectees = 0, nInfectors = 0, Ri = 0))
  )
  
  
  #### Creating day, month, and year variables ####
  
  #Finding all days, months, and years represented in the outbreak
  days <- seq(min(repNumI$observationDate), max(repNumI$observationDate), by = "days")
  daysDf <- cbind.data.frame(date = format(days, "%Y-%m-%d"),
                             dayR = 1:length(days), stringsAsFactors = FALSE)
  
  months <- format(seq(min(repNumI$observationDate), max(repNumI$observationDate), by = "months"), "%Y-%m")
  monthsDf <- cbind.data.frame(month = months, monthR = 1:length(months), stringsAsFactors = FALSE)
  
  years <- year(seq(min(repNumI$observationDate), max(repNumI$observationDate), by = "years"))
  yearsDf <- cbind.data.frame(year = years, yearR = 1:length(years), stringsAsFactors = FALSE)

  
  #Combining the year and month rankings with the raw data
  repNumI2 <- (repNumI
               %>% mutate(month = format(observationDate, "%Y-%m"),
                          year = year(observationDate),
                          date = format(observationDate, "%Y-%m-%d"))
               %>% left_join(daysDf, by = "date")
               %>% left_join(monthsDf, by = "month")
               %>% left_join(yearsDf, by = "year")
               #Making sure monthR starts at 1
               #%>% mutate(monthR = monthR - min(monthR) + 1)
  )
  
  
  
  return(repNumI2)
}


