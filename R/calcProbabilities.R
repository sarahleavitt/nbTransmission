
#' Calculates Relative Transmission Probabilities
#'
#' The function \code{nbProbabilities} uses naive Bayes and an interative estimation
#' procedure to calculate relative transmission probabilities
#'
#' This algorithm takes a dataset of ordered possible infector-infectee pairs in an
#' infectious disease cluster and estimates the relative probability the cases are
#' linked by direct transmission using a classification technique called naive Bayes (NB).
#' 
#' NB is a simple machine learning algorithm that uses Bayes rule to estimate the
#' probability of an outcome in a prediction dataset given a set of covariates from
#' the observed frequencies in a training dataset.
#' 
#' A subset of cases with pathogen WGS or contact investigation data are used to create
#' a training dataset of probable links and non/links. These probable links and 
#' non/links are defined by \code{goldStdVar} which should be a logical vector with
#' \code{TRUE} indicating links, \code{FALSE} nonlinks, and \code{NA} if
#' the pair cannot be used to train (does not have the information or is indeterminate).
#' The covariates can be any categorical variables and could represent
#' spatial, clinical, demographic, and temporal characteristics of the cases. 
#' 
#' Because the outcomes in the training set represent probable and not certain 
#' transmission events and a given case could have mulitple probable infectors, 
#' we use an iterative estimation procedure. This procedure randomly chooses one
#' link of all of the possible links to include in the training dataset \code{nReps}
#' times, and then uses \code{mxn} cross prediction to give all pairs a turn 
#' in the prediction dataset.
#'
#' @param orderedPair The name of the ordered pair-level dataset with the covariates.
#' @param indIDVar The variable name (in quotes) of the individual ID variables 
#' (dataframe \code{orderedPair} must have variables called \code{<indIDVar>.1} and \code{<indIDVar>.2}).
#' @param pairIDVar The variable name (in quotes) of the pair ID variable.
#' @param goldStdVar The variable name (in quotes) that of logical vector defining training links/non-links
#' @param covariates A character vector containing the covariate variable names (in quotes).
#' All covariates need to be categorical factor variables.
#' @param label An optional label string for the run.
#' @param l Laplace smoothing parameter that is added to each cell.
#' @param n The number of folds for nxm cross prediction.
#' @param m The number of times to create n folds in nxm cross prediction.
#' @param nReps The number of times to randomly select the "true" infector.
#'
#' @return List containing two dataframes:
#' \enumerate{
#'   \item \code{probabilities} - a dataframe of transmission probabilities. Column names:
#'      \itemize{
#'        \item \code{label} - the optional label of the run.
#'        \item \code{<pairIDVar>} - the pair ID with the name specified.
#'        \item \code{pAvg} - the mean transmission probability for the pair over all runs.
#'        \item \code{pSD} - the standard deviation of the transmission probability for the pair
#'         over all runs.
#'        \item \code{pScaled} - the mean relative transmission probability for the pair over.
#'         all runs: pAvg scaled so that the probabilities for all infectors per infectee add to 1.
#'        \item \code{pRank} - the rank of the probability of the the pair out of all pairs for that
#'        infectee (in case of ties all values have the minimum rank of the group).
#'        \item \code{nSamples} - the number of probability estimates that contributed to pAvg. This
#'        represents the number of prediction datasets this pair was included in over the \code{nxm}
#'        cross prediction repeated \code{nReps} times.
#'      }
#'   \item \code{estimates} - a dataframe with the effect estimates. Column names:
#'      \itemize{
#'        \item \code{level} - the covariate name and level
#'        \item \code{label} - the optional label of the run
#'        \item \code{oddsMean} - the mean value of the likelihood odds across runs
#'        \item \code{oddsMin} - the min value of the likelihood odds across runs
#'        \item \code{oddsMax} - the max value of the likelihood odds across runs
#'        \item \code{oddsSD} - the standard deviation of the likelihood odds across runs
#'        \item \code{nSamples} - the number of samples included in the average: \code{n*m*nReps}
#'      }
#' }
#'
#' @examples
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
#' #NOTE should run with nReps > 1.
#' covariates = c("Z1", "Z2", "Z3", "Z4", "timeCat")
#' resGen <- nbProbabilities(orderedPair = orderedPair,
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
#' @export



nbProbabilities <- function(orderedPair, indIDVar, pairIDVar, goldStdVar, covariates,
                            label = "", l = 1, n = 10, m = 1, nReps = 10){
  
  orderedPair <- as.data.frame(orderedPair)
  #Creating variables with the individual indID variable
  indIDVar1 <- paste0(indIDVar, ".1")
  indIDVar2 <- paste0(indIDVar, ".2")
  
  #Checking that the named variables are in the dataframe
  if(!indIDVar1 %in% names(orderedPair)){
    stop(paste0(indIDVar1, " is not in the dataframe."))
  }
  if(!indIDVar2 %in% names(orderedPair)){
    stop(paste0(indIDVar2, " is not in the dataframe."))
  }
  if(!pairIDVar %in% names(orderedPair)){
    stop(paste0(pairIDVar, " is not in the dataframe."))
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
  
  #Subsetting to only relevant columns to be more efficient
  orderedPair <- orderedPair[, c(indIDVar1, indIDVar2, pairIDVar,
                                 goldStdVar, covariates)]
  
  #Finding all pairs that can be included in the training dataset
  #And subsetting into the potential links
  posTrain <- orderedPair[!is.na(orderedPair[, goldStdVar]), ]
  posLinks <- posTrain[posTrain[, goldStdVar] == TRUE, ]
  
  

  #### Cross-prediction Procedure ####
  
  #Initializing dataframes to hold results and coefficients
  rAll <- data.frame("p" = numeric(), pairIDVar = character())
  names(rAll) <- c("p", pairIDVar)
  cAll <- data.frame("level" = character(), "odds" = numeric())

  
  pb <- utils::txtProgressBar(min = 0, max = nReps, style = 3)
  for (k in 1:nReps){
    
    #Randomly choosing the "true" infector from all possible
    #Calculating probabilities using mxn cross prediction
    cvResults <- runCV(posTrain, posLinks, orderedPair,
                       indIDVar, pairIDVar, goldStdVar, 
                       covariates, l, n, m)
    rAll <- rbind(rAll, cvResults$rFolds)
    cAll <- rbind(cAll, cvResults$cFolds)
    utils::setTxtProgressBar(pb, k)
  }
  
  
  #### Summarizing Probabilities Over Runs ####
  
  #Averaging the probabilities over all the replicates
  sumData1 <- dplyr::group_by(rAll, !!rlang::sym(pairIDVar))
  sumData2 <- dplyr::summarize(sumData1,
                               pAvg = mean(!!rlang::sym("p"), na.rm = TRUE),
                               pSD = stats::sd(!!rlang::sym("p"), na.rm = TRUE),
                               nSamples = sum(!is.na(!!rlang::sym("p"))),
                               label = dplyr::first(label))
  sumData2 <- dplyr::ungroup(sumData2)
  
  probs <- as.data.frame(dplyr::full_join(sumData2, orderedPair, by = pairIDVar))
  
  #Calculating the total of all probabilities per infectee
  totalP <- stats::aggregate(probs$pAvg, by = list(probs[, indIDVar2]), sum, na.rm = TRUE)
  names(totalP) <- c(indIDVar2, "pTotal")

  #Calculating the scaled probabilities
  probs2 <- merge(probs, totalP, by = indIDVar2)
  probs2$pScaled <- ifelse(probs2$pTotal != 0, probs2$pAvg / probs2$pTotal, 0)
  #Ranking the probabilities for each possible infector
  #Ties are set to the minimum rank of that group
  probs2 <- probs2[order(probs2[, indIDVar2], -probs2$pScaled), ]
  probs2$pRank <- stats::ave(-probs2$pScaled, probs2[, indIDVar2], 
                            FUN = function(x){
                              rank(x, ties.method = "min") 
                            })
  
  #Only keeping columns of interest
  probs2 <- probs2[, c("label", pairIDVar, "pAvg", "pSD", "pScaled", "pRank", "nSamples")]
  
  
  #### Summarizing Measures of Effect Over Runs ####
  
  #Averaging over the measures of effect
  coeffL <- by(cAll,
                 INDICES = list(cAll$level),
                 FUN = function(x){
                   data.frame("level" = unique(x$level),
                              "oddsMean" = mean(x$or, na.rm = TRUE),
                              "oddsMin" = min(x$or, na.rm = TRUE),
                              "oddsMax" = max(x$or, na.rm = TRUE),
                              "oddsSD" = stats::sd(x$or, na.rm = TRUE),
                              "nSamples" = sum(!is.na(x$or)),
                              "label" = label)
                 })
  coeff <- do.call(rbind, coeffL)
  
  return(list("probabilities" = probs2, "estimates" = coeff))
}





runCV <- function(posTrain, posLinks, orderedPair, indIDVar, pairIDVar,
                  goldStdVar, covariates, l, n, m){
  
  #Creating variables with the individual indID variable
  indIDVar1 <- paste0(indIDVar, ".1")
  indIDVar2 <- paste0(indIDVar, ".2")
  
  #Choosing the true infector from all possibles (if multiple)
  #Then subsetting to complete pairs, grouping by infectee, and randomly choosing
  #one possible infector
  linksL <- by(posLinks,
               INDICES = list(posLinks[, indIDVar2]),
               FUN = function(x){
                 x[sample(nrow(x), 1), ]
               })
  links <- do.call(rbind, linksL)
  links$linked <- TRUE
  links <- links[, c(pairIDVar, indIDVar2, "linked")]
  
  #Combining the links with the non-links that do not share an infectee with the links
  trainingFull <- merge(orderedPair, links, by = c(pairIDVar, indIDVar2), all = TRUE)
  trainingFull2 <- trainingFull[trainingFull[, pairIDVar] %in% links[, pairIDVar] |
                                  (trainingFull[, pairIDVar] %in% posTrain[, pairIDVar] &
                                     !trainingFull[, indIDVar2] %in% links[, indIDVar2]), ]
  trainingFull2[is.na(trainingFull2$linked), "linked"] <- FALSE
  

  
  #Creating the cross-valindIDation folds for that part of the training dataset
  cv_splits <- caret::createMultiFolds(trainingFull2$linked, k = n, times = m)
  
  #Initializing dataframes to hold results and coefficients
  rFolds <- data.frame("p" = numeric(), pairIDVar = character())
  names(rFolds) <- c("p", pairIDVar)
  cFolds <- data.frame("level" = character(), "odds" = numeric())
  
  #Running the methods for all of the CV Folds
  for (i in 1:length(cv_splits)){
    
    #Finding training dataset
    trainingPairID <- trainingFull2[, pairIDVar][cv_splits[[i]]]
    trainingRaw <- trainingFull2[trainingFull2[, pairIDVar] %in% trainingPairID, ]
    
    #Finding the infectee whose infector was found in the training dataset
    foundInfector <- trainingRaw[trainingRaw$linked == TRUE, ]
    #Finding all pairs that share an infectee with the pairs where the infector was found
    shareInfectee <- orderedPair[!orderedPair[, pairIDVar] %in% foundInfector[, pairIDVar] &
                                   orderedPair[, indIDVar2] %in% foundInfector[, indIDVar2], ]
    shareInfectee$linked <- FALSE

    #Creating the training datasets
    #Setting probabilities to 1 for training links and 0 for training non-links that are
    #defined as such by the gold standard (not just by sharing an infectee in the above df)
    training <- rbind(trainingRaw, shareInfectee)
    training$p <- ifelse(training[, goldStdVar] == FALSE, 0,
                  ifelse(training[, "linked"] == TRUE, 1, NA))

    #Creating the prediction dataset
    prediction <- dplyr::full_join(orderedPair, links, by = c(pairIDVar, indIDVar2))
    prediction <- prediction[!prediction[, pairIDVar] %in% training[, pairIDVar], ]
    prediction[is.na(prediction$linked), "linked"] <- FALSE
    
    #Calculating probabilities for one split
    sim <- performNB(training, prediction, obsIDVar = pairIDVar,
                     goldStdVar = "linked", covariates, l)
    
    #Combining the results from fold run with the previous folds
    rFolds <- rbind(rFolds, sim[[1]])
    cFolds <- rbind(cFolds, sim[[2]])
  }
  
  return(list("rFolds" = rFolds, "cFolds" = cFolds))
}
