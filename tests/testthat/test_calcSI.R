
context("estimateSI Function")
library(nbTransmission)

#Creating a function with defaults equal to my simulated data
estimateSIWrapper <- function(nbResults,
                             indIDVar = "individualID",
                             timeDiffVar = "infectionDiffY",
                             pVar = "pScaled",
                             clustMethod = c("none", "n", "kd", "hc_absolute", "hc_relative"),
                             cutoffs = NULL,
                             initialPars = c(2,2),
                             bootSamples = 0,
                             alpha = 0.05,
                             epsilon = 0.0001){
  
  siData <- estimateSI(nbResults, indIDVar = indIDVar,
                       timeDiffVar = timeDiffVar, pVar = pVar, clustMethod = clustMethod,
                       cutoffs = cutoffs, initialPars = initialPars, epsilon = epsilon,
                       bootSamples = bootSamples, alpha = alpha)
  
  return(siData)
}

#Shortening dataset for sake of speed
testData <- nbResults[1:500, ]

#Run with no clustering and no CI
siNone <- estimateSIWrapper(testData, clustMethod = "none")
#Run with kd clustering and no CI
siKD <- estimateSIWrapper(testData, clustMethod = "kd", cutoffs = c(0.01, 0.02))
#Run with hc_absolute clustering and no CI
siHC <- estimateSIWrapper(testData, clustMethod = "hc_absolute", cutoffs = c(0.02, 0.05))
#Run with hc_relative clustering and no CI
siHCR <- estimateSIWrapper(testData, clustMethod = "hc_relative", cutoffs = c(2, 3))
#Run with n clustering and no CI
siN <- estimateSIWrapper(testData, clustMethod = "n", cutoffs = 1)

#Run with no clustering and CI
siNoneCI <- estimateSIWrapper(testData, clustMethod = "none", bootSamples = 2)
#Run with clustering and CI
siHCCI <- estimateSIWrapper(testData, clustMethod = "hc_absolute", cutoff = 0.05,
                            bootSamples = 2)



test_that("estimateSI returns a dataframe with the right number of rows",{
  
  expect_true(is.data.frame(siNone))
  expect_true(is.data.frame(siKD))
  expect_true(is.data.frame(siHC))
  expect_true(is.data.frame(siHCR))
  expect_true(is.data.frame(siN))
  expect_true(is.data.frame(siNoneCI))
  expect_true(is.data.frame(siHCCI))
  
  expect_true(nrow(siNone) == 1)
  expect_true(nrow(siKD) == 3)
  expect_true(nrow(siHC) == 3)
  expect_true(nrow(siHCR) == 3)
  expect_true(nrow(siN) == 1)
  expect_true(nrow(siNoneCI) == 1)
  expect_true(nrow(siHCCI) == 1)
  
})

test_that("estimateSI with returns the right column names based on bootSamples",{
  
  expect_true("meanCILB" %in% names(siNoneCI))
  expect_true("meanCILB" %in% names(siHCCI))
  
  expect_false("meanCILB" %in% names(siNone))
  expect_false("meanCILB" %in% names(siKD))
  expect_false("meanCILB" %in% names(siHC))
})



test_that("Descriptive error messages returned",{
  
  expect_error(estimateSIWrapper(nbResults, indIDVar = "garbage"),
               "garbage.1 is not in the data frame.")
  
  expect_error(estimateSIWrapper(nbResults, timeDiffVar = "garbage"),
               "garbage is not in the data frame.")
  
  expect_error(estimateSIWrapper(nbResults, pVar = "garbage"),
               "garbage is not in the data frame.")
  
  #Removing individualID columns
  nbResults2 <- nbResults[!names(nbResults) %in% c("individualID.1")]
  expect_error(estimateSIWrapper(nbResults2, indIDVar = "individualID"),
               "individualID.1 is not in the data frame.")
  
  nbResults3 <- nbResults[!names(nbResults) %in% c("individualID.2")]
  expect_error(estimateSIWrapper(nbResults3, indIDVar = "individualID"),
               "individualID.2 is not in the data frame.")
  
  #Testing that missing clustMethod gets set to none and gives warning
  expect_warning(siMiss <- estimateSIWrapper(testData),
                 "No clustMethod was provided so it was set to 'none'")
  expect_true(siMiss$clustMethod == "none")
  
  #Providing an invalid clustering method
  expect_error(estimateSIWrapper(nbResults, clustMethod = "garbage"),
               "clustMethod must be one of: none, n, kd, hc_absolute, hc_relative")
  
  #Providing a clust method with no cutoff
  expect_error(estimateSIWrapper(nbResults, clustMethod = "hc_absolute"),
               "Please provide one or more cutoff values")
  
  #Too few infectors
  expect_message(estimateSIWrapper(testData, clustMethod = "hc_absolute", cutoff = 0.5),
                 "hc_absolute and 0.5, fewer than 10 individuals would be used for estimation")
  
})


