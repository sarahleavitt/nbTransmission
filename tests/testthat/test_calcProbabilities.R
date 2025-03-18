library(testthat)

context("new nbProbabilities Function")
library(nbTransmission)

## Use the pairData dataset which represents a TB-like outbreak
# First create a dataset of ordered pairs
# orderedPair <- readRDS("covarOrderedPair.rds")
orderedPair <- pairData[pairData$infectionDiffY >= 0, ]

## Create a variable called snpClose that will define probable links
# (<3 SNPs) and nonlinks (>12 SNPs) all pairs with between 2-12 SNPs
# will not be used to train.
orderedPair$snpClose <- ifelse(orderedPair$snpDist < 3, TRUE,
                               ifelse(orderedPair$snpDist > 12, FALSE, NA))

## Running the algorithm

#Creating a function with defaults equal to my simulated data
nbProbWrapper <- function(orderedPair,
                          indIDVar = "individualID",
                          pairIDVar = "pairID",
                          goldStdVar = "snpClose",
                          covariates = c("Z1", "Z2", "Z3", "Z4", "timeCat"),
                          nReps = 1,
                          orType = "adjusted",
                          nBS = 5,
                          pSampled = 1,
                          progressBar = FALSE){
  resGen <- nbProbabilities(orderedPair = orderedPair,
                               indIDVar = indIDVar,
                               pairIDVar = pairIDVar,
                               goldStdVar = goldStdVar,
                               covariates = covariates,
                               nReps = nReps,
                               orType = orType,
                               nBS = nBS,
                               pSampled = pSampled,
                               progressBar = progressBar)
  return(resGen)
}


test_that("nbProbabilities returns a list of two data frames for valid input",{

  resGen <- nbProbWrapper(orderedPair)
  expect_true(is.data.frame(resGen[[1]]))
  expect_true(is.data.frame(resGen[[2]]))

  resGen2 <- nbProbWrapper(orderedPair, progressBar = TRUE)
  expect_true(is.data.frame(resGen2[[1]]))
  expect_true(is.data.frame(resGen2[[2]]))
})


test_that("Descriptive error messages returned",{

  expect_error(nbProbWrapper(orderedPair, indIDVar = "garbage"),
               "garbage.1 is not in the data frame.")

  expect_error(nbProbWrapper(orderedPair, pairIDVar = "garbage"),
               "garbage is not in the data frame.")

  expect_error(nbProbWrapper(orderedPair, goldStdVar = "garbage"),
               "garbage is not in the data frame.")

  expect_error(nbProbWrapper(orderedPair, orType = "garbage"),
               "orType must be either 'adjusted' or 'univariate'.")

  #Removing individualID.1 column
  orderedPair2 <- orderedPair[!names(orderedPair) %in% c("individualID.1")]
  expect_error(nbProbWrapper(orderedPair2, indIDVar = "individualID"),
               "individualID.1 is not in the data frame.")

  #Removing individualID.1 column
  orderedPair3 <- orderedPair[!names(orderedPair) %in% c("individualID.2")]
  expect_error(nbProbWrapper(orderedPair3, indIDVar = "individualID"),
               "individualID.2 is not in the data frame.")

  #Using garbage covariates
  expect_error(nbProbWrapper(orderedPair, covariates = "garbage"),
               "At least one of the covariates is not in the data frame.")

  #Making some covariates factors
  orderedPair4 <- orderedPair
  orderedPair4$Z1 <- as.character(orderedPair4$Z1)
  expect_error(nbProbWrapper(orderedPair4),
               paste0("Z1 is not a factor"))

  orderedPair4$Z2 <- as.character(orderedPair4$Z2)
  expect_error(nbProbWrapper(orderedPair4),
               paste0("Z1, Z2 are not factors"))

})


