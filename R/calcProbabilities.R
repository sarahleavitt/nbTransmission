
#' Calculates Relative Transmission Probabilities
#'
#' Uses naive Bayes and cross vaindIDation to calculate the relative transmission probabilities
#'
#' Add details section
#'
#' @param orderedPair The name of the ordered pair-level dataset with the covariates.
#' @param indIDVar The variable name (in quotes) of the individual ID variable.
#' @param edgeIDVar The variable name (in quotes) of the edge ID variable.
#' @param goldStdVar The variable name (in quotes) that will define linking status.
#' @param covariates A character vector containing the covariate variable names.
#' @param label An optional label string for the run (default is NULL).
#' @param nbWeighting A logical scalar. Do you want to use deep frequency weighting in NB (default is FALSE)?
#' @param n The number of folds for nxm cross valindIDation.
#' @param m The number of times to create n folds in nxm cross valindIDation.
#' @param nReps The number of times to randomly select one infector.
#'
#' @return List containing two dataframes: "probs" with pairdata with an extra column with the
#'    average probabilities over all of the cross valindIDation runs and "coeff" with the
#'    average, minimum, and maximum coefficient values over all cv runs.
#'
#' @examples
#' #Insert example here
#'
#' @import dplyr
#' 
#' @export


calcProbabilities <- function(orderedPair, indIDVar, edgeIDVar, goldStdVar,
                              covariates, label = "", nbWeighting = FALSE,
                              n = 10, m = 1, nReps = 10){
  
  #### Setting up data frames ####
  
  #Creating correctly named variables
  orderedPair$indID.1 <- orderedPair[, paste0(indIDVar, ".1")]
  orderedPair$indID.2 <- orderedPair[, paste0(indIDVar, ".2")]
  orderedPair$edgeID <- orderedPair[, edgeIDVar]
  orderedPair$goldStd <- orderedPair[, goldStdVar]
  
  #Finding all pairs that can be included in the training dataset
  #And subsetting into the potential links
  posTrain <- orderedPair %>% filter(!is.na(goldStd))
  posLinks <- posTrain %>% filter(goldStd == TRUE)
  

  #### Cross-ValindIDation Procedure ####
  
  #Initializing dataframes to hold results and coefficients
  rAll <- NULL
  cAll <- NULL

  for (k in 1:nReps){
    
    #Randomly choosing the "true" infector from all possible
    #Calculating probabilities using mxn cross validation
    cvResults <- runCV(posTrain, posLinks, orderedPair, covariates, goldStdVar,
                       nbWeighting, label, n, m)
    rAll <- bind_rows(rAll, cvResults$rFolds)
    cAll <- bind_rows(cAll, cvResults$cFolds)
  }
  
  
  #### Summarizing Over Runs ####
  
  #Averaging the probabilities over all the replicates
  probs <- (rAll
            %>% group_by(edgeID)
            %>% summarize(pAvg = mean(p, na.rm = TRUE),
                          pSD = sd(p, na.rm = TRUE),
                          nSamples = sum(!is.na(p)),
                          label = first(label))
            %>% full_join(orderedPair, by = "edgeID")
            %>% ungroup()
  )
  
  #Calculating scaled probabilities
  totalP <- (probs
             %>% group_by(indID.2)
             %>% summarize(pTotal = sum(pAvg, na.rm = TRUE))
  )
  probs2 <- (probs
             %>% full_join(totalP, by = "indID.2")
             %>% mutate(pScaled = ifelse(pTotal != 0, pAvg / pTotal, 0))
             %>% select(label, edgeID, pAvg, pSD, pScaled, nSamples)
  )
  
  #Averaging over the measures of effect
  coeff <- (cAll
            %>% group_by(level)
            %>% summarize(label = first(label),
                          ratioMean = mean(ratio, na.rm = TRUE),
                          ratioMin = min(ratio, na.rm = TRUE),
                          ratioMax = max(ratio, na.rm = TRUE),
                          ratioSD = sd(ratio, na.rm = TRUE))
            %>% ungroup()
  )
  
  return(list("probabilities" = probs2, "estimates" = coeff))
}





runCV <- function(posTrain, posLinks, orderedPair, covariates, goldStdVar, nbWeighting, label, n, m){
  
  #Choosing the true infector from all possibles (if multiple)
  #Then subsetting to complete pairs, grouping by infectee, and randomly choosing
  #one possible infector
  links <- (posLinks
            %>% group_by(indID.2)
            %>% sample_n(1)
            %>% ungroup(indID.2)
            %>% mutate(linked = TRUE)
            %>% select(edgeID, indID.2, linked)
  )
  #Combining the links with the non-links that do not share an infectee with the links
  trainingFull <- (orderedPair
                   %>% full_join(links, by = c("edgeID", "indID.2"))
                   %>% filter(edgeID %in% links$edgeID | 
                                (edgeID %in% posTrain$edgeID &
                                   !indID.2 %in% links$indID.2))
                   %>% tidyr::replace_na(list(linked = FALSE))
  )
  
  #Creating the cross-valindIDation folds for that part of the training dataset
  cv_splits <- caret::createMultiFolds(trainingFull$linked, k = n, times = m)
  
  #Initializing dataframes to hold results and coefficients
  rFolds <- NULL
  cFolds <- NULL
  
  #Running the methods for all of the CV Folds
  for (i in 1:length(cv_splits)){
    
    #Finding training dataset
    trainingPairID <- trainingFull$edgeID[cv_splits[[i]]]
    trainingRaw <- trainingFull %>% filter(edgeID %in% trainingPairID)
    #Finding the infectee whose infector was found in the training dataset
    foundInfector <- trainingRaw[trainingRaw$linked == TRUE, ]
    #Finding all pairs that share an infectee with the pairs where the infector was found
    shareInfectee <- (orderedPair
                      %>% filter(!edgeID %in% foundInfector$edgeID &
                                   indID.2 %in% foundInfector$indID.2)
                      %>% mutate(linked = FALSE)
    )
    training <- (trainingRaw
                 %>% bind_rows(shareInfectee)
                 %>% mutate(p = ifelse(goldStd == FALSE, 0,
                                ifelse(linked == TRUE, 1, NA)))
    )
    validation <- (orderedPair
                      %>% full_join(links, by = c("edgeID", "indID.2"))
                      %>% filter(!edgeID %in% training$edgeID)
                      %>% tidyr::replace_na(list(linked = FALSE))
    )
    
    #Calculating probabilities for one split
    sim <- performNB(training, validation, covariates,
                     goldStdVar, nbWeighting, label)
    
    #Combining the results from fold run with the previous folds
    rFolds <- bind_rows(rFolds, sim[[1]])
    cFolds <- bind_rows(cFolds, sim[[2]])
  }
  
  return(list("rFolds" = rFolds, "cFolds" = cFolds))
}
