
context("clusterInfectors Function")
library(nbTransmission)


#Creating a function with defaults equal to my simulated data
clusterInfWrapper <- function(nbResults,
                             indIDVar = "individualID",
                             pVar = "pScaled",
                             clustMethod = c("n", "kd", "hc_absolute", "hc_relative"),
                             cutoff = NA){
  
  siData <- clusterInfectors(nbResults, indIDVar = indIDVar, pVar = pVar,
                             clustMethod = clustMethod, cutoff = cutoff)
  
  return(siData)
}


clustKD1 <- clusterInfWrapper(nbResults, clustMethod = "kd", cutoff = 0.02)
clustKD2 <- clusterInfWrapper(nbResults, clustMethod = "kd", cutoff = 0.1)
clustHC1 <- clusterInfWrapper(nbResults, clustMethod = "hc_absolute", cutoff = 0.01)
clustHC2 <- clusterInfWrapper(nbResults, clustMethod = "hc_absolute", cutoff = 0.1)
clustHCR <- clusterInfWrapper(nbResults, clustMethod = "hc_relative", cutoff = 2)
clustN <- clusterInfWrapper(nbResults, clustMethod = "n", cutoff = 1)



test_that("clusterInf returns a dataframe with cluster column",{
  
  expect_true(nrow(clustKD1) == nrow(nbResults))
  expect_true(nrow(clustKD2) == nrow(nbResults))
  expect_true(nrow(clustHC1) == nrow(nbResults))
  expect_true(nrow(clustHC2) == nrow(nbResults))
  expect_true(nrow(clustHCR) == nrow(nbResults))
  expect_true(nrow(clustN) == nrow(nbResults))
  
  expect_true("cluster" %in% names(clustKD1))
  expect_true("cluster" %in% names(clustKD2))
  expect_true("cluster" %in% names(clustHC1))
  expect_true("cluster" %in% names(clustHC2))
  expect_true("cluster" %in% names(clustHCR))
  expect_true("cluster" %in% names(clustN))
  
})



test_that("Descriptive error messages returned",{
  
  expect_error(clusterInfWrapper(nbResults, indIDVar = "garbage"),
               "garbage.1 is not in the data frame.")
  
  expect_error(clusterInfWrapper(nbResults, pVar = "garbage"),
               "garbage is not in the data frame.")
  
  #Removing individualID columns
  nbResults2 <- nbResults[!names(nbResults) %in% c("individualID.1")]
  expect_error(clusterInfWrapper(nbResults2, indIDVar = "individualID"),
               "individualID.1 is not in the data frame.")
  
  nbResults3 <- nbResults[!names(nbResults) %in% c("individualID.2")]
  expect_error(clusterInfWrapper(nbResults3, indIDVar = "individualID"),
               "individualID.2 is not in the data frame.")
  
  #Testing that missing clustMethod gives an error
  expect_error(clusterInfWrapper(nbResults),
                 "Please provide a clustering method")
  
  #Providing an invalid clustering method
  expect_error(clusterInfWrapper(nbResults, clustMethod = "garbage"),
               "clustMethod must be one of: n, kd, hc_absolute, hc_relative")
  
})
