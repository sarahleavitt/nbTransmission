
context("nbProbabilities Function")
library(nbTransmission)

## Use the pairData dataset which represents a TB-like outbreak
# First create a dataset of ordered pairs
orderedPair <- pairData[pairData$infectionDiffY > 0, ]

## Create a variable called snpClose that will define probable links
# (<3 SNPs) and nonlinks (>12 SNPs) all pairs with between 2-12 SNPs
# will not be used to train.
orderedPair$snpClose <- ifelse(orderedPair$snpDist < 3, TRUE,
                               ifelse(orderedPair$snpDist > 12, FALSE, NA))
covariates = c("Z1", "Z2", "Z3", "Z4", "timeCat")

## Running the algorithm

#Creating a function with defaults equal to my simulated data
nbProbWrapper <- function(orderedPair,
                            indIDVar = "individualID",
                            pairIDVar = "pairID",
                            goldStdVar = "snpClose",
                            covariates = c("Z1", "Z2", "Z3", "Z4", "timeCat"),
                            nReps = 1){
  
  resGen <- nbProbabilities(orderedPair = orderedPair,
                              indIDVar = indIDVar,
                              pairIDVar = pairIDVar,
                              goldStdVar = goldStdVar,
                              covariates = c("Z1", "Z2", "Z3", "Z4", "timeCat"),
                              nReps = nReps)
  return(resGen)
}


test_that("nbProbabilities returns a list of two data frames for valid input",{
  
  resGen <- nbProbWrapper(orderedPair)
  expect_true(is.data.frame(resGen[[1]]))
  expect_true(is.data.frame(resGen[[2]]))
})


test_that("Descriptive error messages returned",{

 expect_error(nbProbWrapper(orderedPair, indIDVar = "garbage"),
              "garbage.1 is not in the data frame.")
 
 #Removing individualID.2 column
 orderedPair2 <- orderedPair[!names(orderedPair) %in% c("individualID.2")]
 expect_error(nbProbWrapper(orderedPair2, indIDVar = "individualID"),
              "individualID.2 is not in the data frame.")
 
 expect_error(nbProbWrapper(orderedPair, pairIDVar = "garbage"),
              "garbage is not in the data frame.")
 
 expect_error(nbProbWrapper(orderedPair, goldStdVar = "garbage"),
              "garbage is not in the data frame.")
  
})

