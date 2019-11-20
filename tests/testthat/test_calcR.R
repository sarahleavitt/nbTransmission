
context("calcR Suite of Functions")
library(nbTransmission)

## Use the pairData dataset which represents a TB-like outbreak
# First create a dataset of ordered pairs
orderedPair <- pairData[pairData$infectionDiffY > 0, ]

## Create a variable called snpClose that will define probable links
# (<3 SNPs) and nonlinks (>12 SNPs) all pairs with between 2-12 SNPs
# will not be used to train.
orderedPair$snpClose <- ifelse(orderedPair$snpDist < 3, TRUE,
                               ifelse(orderedPair$snpDist > 12, FALSE, NA))
table(orderedPair$snpClose)

## Running the algorithm
covariates = c("Z1", "Z2", "Z3", "Z4", "timeCat")
resGen <- nbProbabilities(orderedPair = orderedPair, indIDVar = "individualID",
                            pairIDVar = "pairID", goldStdVar = "snpClose",
                            covariates = covariates, nReps = 1)
allProbs <- merge(resGen[[1]], orderedPair, by = "pairID", all = TRUE)


#Creating a function with defaults equal to my simulated data
estimateRWrapper <- function(allProbs,
                         dateVar = "infectionDate",
                         indIDVar = "individualID",
                         pVar = "pScaled",
                         timeFrame = "months",
                         rangeForAvg = NULL,
                         bootSamples = 0,
                         alpha = 0.05){
  
  rList <- estimateR(allProbs, dateVar = dateVar, indIDVar = indIDVar,
              pVar = pVar, timeFrame = timeFrame, rangeForAvg = rangeForAvg,
              bootSamples = bootSamples, alpha = alpha)
  
  return(rList)
}

#Run with range specified and no CI
rDataR <- estimateRWrapper(allProbs, rangeForAvg = c(10, 100))
#Run with range specified and CI
rDataRC <- estimateRWrapper(allProbs, rangeForAvg = c(10, 100), bootSamples = 2)

#Run with no range specified and no CI
#rData <- estimateRWrapper(allProbs)
#Run with no range specified and CI
#rDataC <- estimateRWrapper(allProbs, bootSamples = 10)


test_that("estimateR returns a list of three data frames for valid input",{
  
  expect_true(is.data.frame(rDataR[[1]]))
  expect_true(is.data.frame(rDataR[[2]]))
  expect_true(is.data.frame(rDataR[[3]]))
  
  expect_true(is.data.frame(rDataRC[[1]]))
  expect_true(is.data.frame(rDataRC[[2]]))
  expect_true(is.data.frame(rDataRC[[3]]))
})

test_that("estimateR with returns the right column names based on bootSamples",{
  
  expect_true("ciLower" %in% names(rDataRC[[2]]))
  expect_true("ciLower" %in% names(rDataRC[[3]]))
  
  expect_false("ciLower" %in% names(rDataR[[2]]))
  expect_false("ciLower" %in% names(rDataR[[3]]))
})


test_that("Descriptive error messages returned",{
  
  expect_error(estimateRWrapper(allProbs, indIDVar = "garbage"),
               "garbage.1 is not in the data frame.")
  
  expect_error(estimateRWrapper(allProbs, dateVar = "garbage"),
               "garbage.1 is not in the data frame.")
  
  expect_error(estimateRWrapper(allProbs, pVar = "garbage"),
               "garbage is not in the data frame.")
  
  #Removing individualID columns
  allProbs2 <- allProbs[!names(allProbs) %in% c("individualID.1")]
  expect_error(estimateRWrapper(allProbs2, indIDVar = "individualID"),
               "individualID.1 is not in the data frame.")

  allProbs3 <- allProbs[!names(allProbs) %in% c("individualID.2")]
  expect_error(estimateRWrapper(allProbs3, indIDVar = "individualID"),
               "individualID.2 is not in the data frame.")
  
  #Removing the date columns
  allProbs4 <- allProbs[!names(allProbs) %in% c("infectionDate.1")]
  expect_error(estimateRWrapper(allProbs4, dateVar = "infectionDate"),
               "infectionDate.1 is not in the data frame.")
  
  allProbs5 <- allProbs[!names(allProbs) %in% c("infectionDate.2")]
  expect_error(estimateRWrapper(allProbs5, dateVar = "infectionDate"),
               "infectionDate.2 is not in the data frame.")
  
  #Changing dates to character variables
  allProbs$infectionDatec.1 <- as.character(allProbs$infectionDate.1)
  allProbs$infectionDatec.2 <- allProbs$infectionDate.2
  expect_error(estimateRWrapper(allProbs, dateVar = "infectionDatec"),
               "infectionDatec.1 must be either a date or a date-time (POSIXt) object.",
               fixed = TRUE)
  
  allProbs$infectionDatec.1 <- allProbs$infectionDate.1
  allProbs$infectionDatec.2 <- as.character(allProbs$infectionDate.2)
  expect_error(estimateRWrapper(allProbs, dateVar = "infectionDatec"),
               "infectionDatec.2 must be either a date or a date-time (POSIXt) object.",
               fixed = TRUE)
  
  allProbs5 <- allProbs[!names(allProbs) %in% c("infectionDate.2")]
  expect_error(estimateRWrapper(allProbs5, indIDVar = "infectionDate"),
               "infectionDate.2 is not in the data frame.")
  
  #Testing timeFrame error
  expect_error(estimateRWrapper(allProbs, timeFrame = "garbage"),
               paste0("timeFrame must be one of: ",
                      paste0( c("days", "months", "weeks", "years"), collapse = ", ")))
})


