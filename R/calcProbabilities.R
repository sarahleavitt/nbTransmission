
#' Calculates Relative Transmission Probabilities
#'
#' \code{calcProbabilities} uses naive Bayes and cross vaindIDation to calculate the relative
#'  transmission probabilities
#'
#' Add details section
#'
#' @param orderedPair The name of the ordered pair-level dataset with the covariates.
#' @param indIDVar The variable name (in quotes) of the individual ID variable.
#' @param edgeIDVar The variable name (in quotes) of the edge ID variable.
#' @param goldStdVar The variable name (in quotes) that will define linking status.
#' @param covariates A character vector containing the covariate variable names.
#' @param label An optional label string for the run (default is NULL).
#' @param l Laplace smoothing parameter that is added to each cell (default is 1).
#' @param nbWeighting A logical scalar. Do you want to use deep frequency weighting in NB (default is FALSE)?
#' @param n The number of folds for nxm cross valindIDation.
#' @param m The number of times to create n folds in nxm cross valindIDation.
#' @param nReps The number of times to randomly select one infector.
#'
#' @return List containing two dataframes:
#' \enumerate{
#'   \item \code{probabilities} - a dataframe of transmission probabilities with the
#'    following columns:
#'      \itemize{
#'        \item \code{label} - the optional label of the run
#'        \item \code{pAvg} - the mean transmission probability for the pair over all runs
#'        \item \code{pSD} - the standard deviation of the transmission probability for the pair
#'         over all runs
#'        \item \code{pScaled} - the mean relative transmission probability for the pair over
#'         all runs: pAvg scaled so that the probabilities for all infectors per infectee add to 1.
#'        \item \code{pRank} - the rank of the probability of the the pair out of all pairs for that
#'        infectee
#'        \item \code{nSamples} - the number of probability estimates that contributed to pAvg. This
#'        represents the number of validation datasets this pair was included in over the \code{nxm}
#'        cross validation repeated \code{nReps} times.
#'        \item \code{edgeIDVar} - the edgeID variable with the name specified
#'      }
#'   \item \code{estimates} - a dataframe with the effect estimates with the following columns:
#'      \itemize{
#'        \item \code{level} - the covariate name and level
#'        \item \code{label} - the optional label of the run
#'        \item \code{ratioMean} - the mean value of the likelihood ratio across runs
#'        \item \code{ratioMin} - the min value of the likelihood ratio across runs
#'        \item \code{ratioMax} - the max value of the likelihood ratio across runs
#'        \item \code{ratioSD} - the standard deviation of the likelihood ratio across runs
#'        \item \code{nSamples} - the number of samples included in the average: \code{n * m * nReps}
#'      }
#' }
#'
#' @examples
#' #Insert example here
#'
#' @import dplyr
#' 
#' @export


calcProbabilities <- function(orderedPair, indIDVar, edgeIDVar, goldStdVar,
                              covariates, label = "", l = 1, nbWeighting = FALSE, 
                              n = 10, m = 1, nReps = 10){
  
  orderedPair <- as.data.frame(orderedPair)
  
  #Checking that the named variables are in the dataframe
  if(!paste0(indIDVar, ".1") %in% names(orderedPair)){
    stop(paste0(paste0(indIDVar, ".1"), " is not in the dataframe."))
  }
  if(!paste0(indIDVar, ".2") %in% names(orderedPair)){
    stop(paste0(paste0(indIDVar, ".2"), " is not in the dataframe."))
  }
  if(!edgeIDVar %in% names(orderedPair)){
    stop(paste0(edgeIDVar, " is not in the dataframe."))
  }
  if(!goldStdVar %in% names(orderedPair)){
    stop(paste0(goldStdVar, " is not in the dataframe."))
  }
  
  #Checking that the covariates are in the dataframe
  covarTest <- covariates %in% names(orderedPair)
  if(FALSE %in% covarTest){
    stop("At least one of the covariates is not in the dataframe.")
  }
  
  #Checking that all of the covariates are factors
  covarDf <- orderedPair[, covariates]
  notFactor <- names(covarDf)[!sapply(covarDf, is.factor)]
  notFactorC <- paste0(notFactor, collapse = ", ")
  if(FALSE %in% sapply(covarDf, is.factor)){
    stop(paste0(notFactorC, " are not a factors"))
  }
  
  
  #### Setting up data frames ####
  
  #Creating correctly named variables
  orderedPair$indID.1 <- orderedPair[, paste0(indIDVar, ".1")]
  orderedPair$indID.2 <- orderedPair[, paste0(indIDVar, ".2")]
  orderedPair$edgeID <- orderedPair[, edgeIDVar]
  orderedPair$goldStd <- orderedPair[, goldStdVar]
  
  #Subsetting to only relevant columns
  orderedPair <- orderedPair[, c("indID.1", "indID.2", "edgeID", "goldStd",
                                 goldStdVar, covariates)]
  
  #Finding all pairs that can be included in the training dataset
  #And subsetting into the potential links
  posTrain <- orderedPair[!is.na(orderedPair$goldStd), ]
  posLinks <- posTrain[posTrain$goldStd == TRUE, ]
  
  

  #### Cross-ValindIDation Procedure ####
  
  #Initializing dataframes to hold results and coefficients
  rAll <- NULL
  cAll <- NULL

  for (k in 1:nReps){
    
    #Randomly choosing the "true" infector from all possible
    #Calculating probabilities using mxn cross validation
    cvResults <- runCV(posTrain, posLinks, orderedPair, covariates, l, nbWeighting, n, m)
    rAll <- bind_rows(rAll, cvResults$rFolds)
    cAll <- bind_rows(cAll, cvResults$cFolds)
  }
  
  
  #### Summarizing Over Runs ####
  
  #Averaging the probabilities over all the replicates
  probs <- (rAll
            %>% group_by(edgeID)
            %>% summarize(pAvg = mean(p, na.rm = TRUE),
                          pSD = stats::sd(p, na.rm = TRUE),
                          nSamples = sum(!is.na(p)))
            %>% mutate(label = label)
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
             #Ranking the probabilities for each possible infector
             #Ties are set to the minimum rank of that group
             %>% group_by(indID.2)
             %>% arrange(desc(pScaled))
             %>% mutate(pRank = rank(desc(pScaled), ties.method = "min"))
             %>% ungroup()
             %>% select(label, edgeID, pAvg, pSD, pScaled, pRank, nSamples)
  )
  #Renaming the edgeID variable to match input
  probs2[, edgeIDVar] <- probs2$edgeID
  if(edgeIDVar != "edgeID"){probs2 <- probs2 %>% select(-edgeID)}
  
  
  #Averaging over the measures of effect
  coeff <- (cAll
            %>% group_by(level)
            %>% summarize(label = first(label),
                          ratioMean = mean(ratio, na.rm = TRUE),
                          ratioMin = min(ratio, na.rm = TRUE),
                          ratioMax = max(ratio, na.rm = TRUE),
                          ratioSD = stats::sd(ratio, na.rm = TRUE),
                          nSamples = sum(!is.na(ratio)))
            %>% mutate(label = label)
            %>% ungroup()
  )
  
  return(list("probabilities" = probs2, "estimates" = coeff))
}





runCV <- function(posTrain, posLinks, orderedPair,
                  covariates, l, nbWeighting, n, m){
  
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
    sim <- performNB(training, validation, edgeIDVar = "edgeID",
                     goldStdVar = "goldStd", covariates, l, nbWeighting)
    
    #Combining the results from fold run with the previous folds
    rFolds <- bind_rows(rFolds, sim[[1]])
    cFolds <- bind_rows(cFolds, sim[[2]])
  }
  
  return(list("rFolds" = rFolds, "cFolds" = cFolds))
}
