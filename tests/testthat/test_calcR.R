
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
resGen <- calcProbabilities(orderedPair = orderedPair, indIDVar = "individualID",
                            pairIDVar = "pairID", goldStdVar = "snpClose",
                            covariates = covariates, nReps = 1)
allProbs <- merge(resGen[[1]], orderedPair, by = "pairID", all = TRUE)


#Creating a function with defaults equal to my simulated data
calcRWrapper <- function(orderedPair,
                         dateVar = "infectionDate",
                         indIDVar = "individualID",
                         pVar = "pScaled",
                         timeFrame = "months",
                         rangeForAvg = NULL,
                         bootSamples = 0,
                         alpha = 0.05){
  
  rList <- calcR(allProbs, dateVar = dateVar, indIDVar = indIDVar,
              pVar = pVar, timeFrame = timeFrame, rangeForAvg = rangeForAvg,
              bootSamples = bootSamples, alpha = alpha)
  
  return(rList)
}

#Run with range specified and no CI
rDataR <- calcRWrapper(allProbs, rangeForAvg = c(10, 100))
#Run with range specified and CI
rDataRC <- calcRWrapper(allProbs, rangeForAvg = c(10, 100), bootSamples = 10)

#Run with no range specified and no CI
#rData <- calcRWrapper(allProbs)
#Run with no range specified and CI
#rDataC <- calcRWrapper(allProbs, bootSamples = 10)


test_that("calcR returns a list of three dataframes for valid input",{
  
  expect_true(is.data.frame(rDataR[[1]]))
  expect_true(is.data.frame(rDataR[[2]]))
  expect_true(is.data.frame(rDataR[[3]]))
  
  expect_true(is.data.frame(rDataRC[[1]]))
  expect_true(is.data.frame(rDataRC[[2]]))
  expect_true(is.data.frame(rDataRC[[3]]))
})

#Test warning for rData and rDataRC
#Test column names contain CIs for RDataC and rDataRC and not for rData and rDataRC


