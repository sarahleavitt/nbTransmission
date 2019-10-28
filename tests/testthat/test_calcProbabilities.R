
context("calcProbabilities Function")
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
calcProbWrapper <- function(orderedPair,
                            indIDVar = "individualID",
                            pairIDVar = "pairID",
                            goldStdVar = "snpClose",
                            covariates = c("Z1", "Z2", "Z3", "Z4", "timeCat"),
                            nReps = 1){
  
  resGen <- calcProbabilities(orderedPair = orderedPair,
                              indIDVar = indIDVar,
                              pairIDVar = pairIDVar,
                              goldStdVar = goldStdVar,
                              covariates = c("Z1", "Z2", "Z3", "Z4", "timeCat"),
                              nReps = nReps)
  return(resGen)
}


test_that("calcProbabilities returns a list of two dataframes for valid input",{
  
  resGen <- calcProbWrapper(orderedPair)
  expect_true(is.data.frame(resGen[[1]]))
  expect_true(is.data.frame(resGen[[2]]))
})


test_that("Descriptive error messages returned",{

 expect_error(calcProbWrapper(orderedPair, indIDVar = "garbage"),
              "garbage.1 is not in the dataframe.")
 
 #Removing individualID.2 column
 orderedPair2 <- orderedPair[!names(orderedPair) %in% c("individualID.2")]
 expect_error(calcProbWrapper(orderedPair2, indIDVar = "individualID"),
              "individualID.2 is not in the dataframe.")
 
 expect_error(calcProbWrapper(orderedPair, pairIDVar = "garbage"),
              "garbage is not in the dataframe.")
 
 expect_error(calcProbWrapper(orderedPair, goldStdVar = "garbage"),
              "garbage is not in the dataframe.")
  
})
