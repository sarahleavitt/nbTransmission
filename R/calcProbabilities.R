
#' Calculate Relative Transmission Probabilities
#'
#' Uses naive Bayes and cross vaidation to calculate the relative transmission probabilities
#'
#' @param indData individual-level dataset with covariate values
#' @param pairData pair-level dataset with covariate values
#' @param observationDate the variable name (in quotes) of the date that the individual is observed
#' @param individualID the variable name (in quotes) of the id variable
#' @param goldStandard the variable name (in quotes) that will define linking status
#' @param covariates vector of the covariate variable names
#' @param scheme optional label for the run (default is NULL)
#' @param n number of folds for nxm cross validation
#' @param m number of times to create n folds in nxm cross validation
#' @param nReps number of times to randomly select one infector
#'
#' @return List containing two dataframes: "probs" with pairdata with an extra column with the probabilities and
#'    "coeff" with the coefficients.
#'
#' @examples
#' #Insert example here
#'
#' @export


calcProbabilities <- function(indData, pairData, observationDate, individualID = "individualID",
                              goldStandard, covariates, scheme = NULL,
                              n = 10, m = 1, nReps = 50){
  
  #Creating correctly named variables
  indData$observationDate <- indData[, observationDate]
  pairData$observationDate.1 <- pairData[, paste0(observationDate, ".1")]
  pairData$observationDate.2 <- pairData[, paste0(observationDate, ".2")]
  pairData$observationDiff <- as.numeric(difftime(pairData$observationDate.2,
                                                  pairData$observationDate.1, units = "days"))
  
  indData$individualID <- indData[, individualID]
  pairData$individualID.1 <- pairData[, paste0(individualID, ".1")]
  pairData$individualID.2 <- pairData[, paste0(individualID, ".2")]
  
  #Creating the edgeID
  pairData <- pairData %>% unite(edgeID, individualID.1, individualID.2, remove = FALSE)
  
  #Restricting to sampled cases
  indData <- indData %>% filter(!is.na(observationDate))
  pairData <- pairData %>% filter(!is.na(observationDiff))
  
  #Subseting to the pairs with the potential infector observed before the infectee
  covarOrderedPair <- pairData %>% filter(observationDiff > 0)
  
  
  
  ############### Finding all possible training pairs ################
  
  #Restricting to cases that can be in the training dataset (no missing variables)
  trainingID <- (indData
                 %>% filter(complete == TRUE)
                 %>% pull(individualID)
  )
  
  #Finding all pairs that involve those people
  posTrain <- covarOrderedPair %>% filter(individualID.1 %in% trainingID &
                                            individualID.2 %in% trainingID)
  #Finding all pairs that can be included in the training dataset (have gold standard)
  posTrain <- posTrain[!is.na(posTrain[, goldStandard]) ,]
  
  #Finding all pairs that can be included in the gold standard dataset as links
  posLinks <- posTrain[posTrain[, goldStandard] == TRUE, ]
  
  
  #### Cross-Validation Probability Calculation ####
  
  #Initializing dataframes to hold results, messages, and performance
  rAll <- NULL
  pAll <- NULL
  cAll <- NULL
  messages <- NULL

  for (k in 1:nReps){
    
    #Choosing the true infector from all possibles (if multiple)
    #Then subsetting to complete pairs, grouping by infectee, and randomly choosing
    #one possible infector
    links <- (posLinks
              %>% group_by(individualID.2)
              %>% sample_n(1)
              %>% ungroup(individualID.2)
              %>% mutate(linked = TRUE)
              %>% select(edgeID, individualID.2, linked)
    )
    #Combining the links with the non-links that do not share an infectee with the links
    trainingFull <- (covarOrderedPair
                     %>% full_join(links, by = c("edgeID", "individualID.2"))
                     %>% filter(edgeID %in% links$edgeID | 
                                  (edgeID %in% posTrain$edgeID &
                                   !individualID.2 %in% links$individualID.2))
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
      shareInfectee <- (covarOrderedPair
                        %>% filter(!edgeID %in% foundInfector$edgeID &
                                     individualID.2 %in% foundInfector$individualID.2)
                        %>% mutate(linked = FALSE)
      )
      training <- trainingRaw %>% bind_rows(shareInfectee)
      validation <- (covarOrderedPair
                     %>% full_join(links, by = c("edgeID", "individualID.2"))
                     %>% filter(!edgeID %in% training$edgeID)
                     %>% replace_na(list(linked = FALSE))
      )
      
      #Calculating probabilities for one split
      # sim1 <- performLogistic(training, validation, covariates, 
      #                        scheme, goldStandard, pTraining = 1, thresholds)
      sim2 <- performBayes(training, validation, covariates, weighting=FALSE,
                           scheme, goldStandard)
      
      #Combining the results from fold run with the previous folds
      rAll <- bind_rows(rAll, sim2[[1]])
      
      #Combining the erros for the various methods
      mTemp <- as.data.frame(sim2[[2]])
      mTemp <- mTemp %>% filter(!is.na(messages)) %>% filter(!duplicated(.))
      messages <- rbind.data.frame(messages, mTemp)
      
      #Extracting significance of coefficients from regular logistic regression
      cAll <- bind_rows(cAll, sim2[[3]])
    }
  }
  
  #Averaging the probabilities over all the replicates
  results <- (rAll
              %>% group_by(label, edgeID)
              %>% summarize(pAvg = mean(p, na.rm = TRUE),
                            pSD = sd(p, na.rm = TRUE),
                            nLinksTrain = mean(nLinksTrain),
                            nSamples = sum(!is.na(p)))
              %>% full_join(covarOrderedPair, by = "edgeID")
  )
  
  #Calculating scaled probabilities
  totalP <- (results
             %>% group_by(label, individualID.2)
             %>% summarize(pTotal = sum(pAvg, na.rm = TRUE))
  )
  probs <- (results
               %>% full_join(totalP, by = c("label", "individualID.2"))
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

