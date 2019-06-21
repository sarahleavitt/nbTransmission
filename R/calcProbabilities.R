
#' Calculates Relative Transmission Probabilities
#'
#' Uses naive Bayes and cross vaidation to calculate the relative transmission probabilities
#'
#' @param indData individual-level dataset with covariate values
#' @param pairData pair-level dataset with covariate values
#' @param dateVar the variable name (in quotes) of the date that the individual is observed
#' @param idVar the variable name (in quotes) of the id variable
#' @param goldStdVar the variable name (in quotes) that will define linking status
#' @param covariates vector of the covariate variable names
#' @param label optional label for the run (default is NULL)
#' @param n number of folds for nxm cross validation
#' @param m number of times to create n folds in nxm cross validation
#' @param nReps number of times to randomly select one infector
#'
#' @return List containing two dataframes: "probs" with pairdata with an extra column with the
#'    average probabilities over all of the cross validation runs and "coeff" with the
#'    average, minimum, and maximum coefficient values over all cv runs.
#'
#' @examples
#' #Insert example here
#'
#' @export


calcProbabilities <- function(indData, pairData, dateVar, idVar,
                              goldStdVar, covariates, label = NULL,
                              n = 10, m = 1, nReps = 50){
  
  #### Setting up data frames ####
  
  #Creating correctly named variables
  indData$date <- indData[, dateVar]
  pairData$date.1 <- pairData[, paste0(dateVar, ".1")]
  pairData$date.2 <- pairData[, paste0(dateVar, ".2")]
  pairData$timeDiff <- as.numeric(difftime(pairData$date.2, pairData$date.1, units = "days"))
  
  indData$id <- indData[, idVar]
  pairData$id.1 <- pairData[, paste0(idVar, ".1")]
  pairData$id.2 <- pairData[, paste0(idVar, ".2")]
  
  #Creating the edgeID
  pairData <- pairData %>% unite(edgeID, id.1, id.2, remove = FALSE)
  
  #Restricting to sampled cases
  indData <- indData %>% filter(!is.na(date))
  pairData <- pairData %>% filter(!is.na(timeDiff))
  
  #Subseting to the pairs with the potential infector observed before the infectee
  orderedPair <- pairData %>% filter(timeDiff > 0)
  
  #Finding all pairs that can be included in the training dataset (have gold standard)
  posTrain <- orderedPair[!is.na(orderedPair[, goldStdVar]) ,]
  
  #Finding all pairs that can be included in the gold standard dataset as links
  posLinks <- posTrain[posTrain[, goldStdVar] == TRUE, ]
  

  
  #### Cross-Validation Procedure ####
  
  #Initializing dataframes to hold results and coefficients
  rAll <- NULL
  cAll <- NULL

  for (k in 1:nReps){
    
    #Choosing the true infector from all possibles (if multiple)
    #Then subsetting to complete pairs, grouping by infectee, and randomly choosing
    #one possible infector
    links <- (posLinks
              %>% group_by(id.2)
              %>% sample_n(1)
              %>% ungroup(id.2)
              %>% mutate(linked = TRUE)
              %>% select(edgeID, id.2, linked)
    )
    #Combining the links with the non-links that do not share an infectee with the links
    trainingFull <- (orderedPair
                     %>% full_join(links, by = c("edgeID", "id.2"))
                     %>% filter(edgeID %in% links$edgeID | 
                                  (edgeID %in% posTrain$edgeID &
                                   !id.2 %in% links$id.2))
                     %>% replace_na(list(linked = FALSE))
    )
    
    #Creating the cross-validation folds for that part of the training dataset
    cv_splits <- createMultiFolds(trainingFull$linked, k = n, times = m)
    
    #Running the methods for all of the CV Folds
    for (i in 1:length(cv_splits)){
      
      #Finding training dataset
      trainingPairID <- trainingFull$edgeID[cv_splits[[i]]]
      trainingRaw <- (trainingFull
                      %>% filter(edgeID %in% trainingPairID)
                      %>% mutate(p = ifelse(linked == TRUE, 1, 0))
      )
      #Finding the infectee whose infector was found in the training dataset
      foundInfector <- trainingRaw[trainingRaw$linked == TRUE, ]
      #Finding all pairs that share an infectee with the pairs where the infector was found
      shareInfectee <- (orderedPair
                        %>% filter(!edgeID %in% foundInfector$edgeID &
                                     id.2 %in% foundInfector$id.2)
                        %>% mutate(linked = FALSE)
      )
      training <- trainingRaw %>% bind_rows(shareInfectee)
      validation <- (orderedPair
                     %>% full_join(links, by = c("edgeID", "id.2"))
                     %>% filter(!edgeID %in% training$edgeID)
                     %>% replace_na(list(linked = FALSE))
      )
      
      #Calculating probabilities for one split
      sim <- performNB(training, validation, covariates, goldStdVar, weighting=FALSE, label)
      
      #Combining the results from fold run with the previous folds
      rAll <- bind_rows(rAll, sim[[1]])
      cAll <- bind_rows(cAll, sim[[2]])
    }
  }
  
  
  #### Summarizing Over Runs ####
  
  #Averaging the probabilities over all the replicates
  results <- (rAll
              %>% group_by(label, edgeID)
              %>% summarize(pAvg = mean(p, na.rm = TRUE),
                            pSD = sd(p, na.rm = TRUE),
                            nLinksTrain = mean(nLinksTrain),
                            nSamples = sum(!is.na(p)))
              %>% full_join(orderedPair, by = "edgeID")
  )
  
  #Calculating scaled probabilities
  totalP <- (results
             %>% group_by(label, id.2)
             %>% summarize(pTotal = sum(pAvg, na.rm = TRUE))
  )
  probs <- (results
               %>% full_join(totalP, by = c("label", "id.2"))
               %>% mutate(pScaled = ifelse(pTotal != 0, pAvg / pTotal, 0))
               %>% select(-pTotal)
               %>% ungroup()
  )
  
  #Averaging over the measures of effect
  coeff <- (cAll
            %>% group_by(label, level)
            %>% summarize(ratioMean = mean(ratio, na.rm = TRUE),
                          ratioMin = min(ratio, na.rm = TRUE),
                          ratioMax = max(ratio, na.rm = TRUE),
                          ratioSD = sd(ratio, na.rm = TRUE))
            %>% ungroup()
  )
  
  return(list(probs, coeff))
}

