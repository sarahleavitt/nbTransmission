
context("estimateR Function")
library(nbTransmission)

#Creating a function with defaults equal to my simulated data
estimateRWrapper <- function(nbResults,
                         dateVar = "infectionDate",
                         indIDVar = "individualID",
                         pVar = "pScaled",
                         timeFrame = "months",
                         rangeForAvg = c(10, 100),
                         bootSamples = 0,
                         alpha = 0.05){
  
  rList <- estimateR(nbResults, dateVar = dateVar, indIDVar = indIDVar,
              pVar = pVar, timeFrame = timeFrame, rangeForAvg = rangeForAvg,
              bootSamples = bootSamples, alpha = alpha)
  
  return(rList)
}

#Run with range specified and no CI
rDataR <- estimateRWrapper(nbResults)
#Run with range specified and CI
rDataRC <- estimateRWrapper(nbResults, bootSamples = 2)


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


test_that("Time frame is correct",{

  expect_true(rDataR$timeFrame == "months")
  
  rDataD <- estimateRWrapper(nbResults, timeFrame = "days")
  expect_true(rDataD$timeFrame == "days")
  
  rDataY <- estimateRWrapper(nbResults, timeFrame = "years")
  expect_true(rDataY$timeFrame == "years")
  
  rDataW <- estimateRWrapper(nbResults, timeFrame = "weeks")
  expect_true(rDataW$timeFrame == "weeks")
  
})


test_that("Descriptive error messages returned",{
  
  expect_error(estimateRWrapper(nbResults, indIDVar = "garbage"),
               "garbage.1 is not in the data frame.")
  
  expect_error(estimateRWrapper(nbResults, dateVar = "garbage"),
               "garbage.1 is not in the data frame.")
  
  expect_error(estimateRWrapper(nbResults, pVar = "garbage"),
               "garbage is not in the data frame.")
  
  #Removing individualID columns
  nbResults2 <- nbResults[!names(nbResults) %in% c("individualID.1")]
  expect_error(estimateRWrapper(nbResults2, indIDVar = "individualID"),
               "individualID.1 is not in the data frame.")

  nbResults3 <- nbResults[!names(nbResults) %in% c("individualID.2")]
  expect_error(estimateRWrapper(nbResults3, indIDVar = "individualID"),
               "individualID.2 is not in the data frame.")
  
  #Removing the date columns
  nbResults4 <- nbResults[!names(nbResults) %in% c("infectionDate.1")]
  expect_error(estimateRWrapper(nbResults4, dateVar = "infectionDate"),
               "infectionDate.1 is not in the data frame.")
  
  nbResults5 <- nbResults[!names(nbResults) %in% c("infectionDate.2")]
  expect_error(estimateRWrapper(nbResults5, dateVar = "infectionDate"),
               "infectionDate.2 is not in the data frame.")
  
  #Changing dates to character variables
  nbResults$infectionDatec.1 <- as.character(nbResults$infectionDate.1)
  nbResults$infectionDatec.2 <- nbResults$infectionDate.2
  expect_error(estimateRWrapper(nbResults, dateVar = "infectionDatec"),
               "infectionDatec.1 must be either a date or a date-time (POSIXt) object.",
               fixed = TRUE)
  
  nbResults$infectionDatec.1 <- nbResults$infectionDate.1
  nbResults$infectionDatec.2 <- as.character(nbResults$infectionDate.2)
  expect_error(estimateRWrapper(nbResults, dateVar = "infectionDatec"),
               "infectionDatec.2 must be either a date or a date-time (POSIXt) object.",
               fixed = TRUE)
  
  nbResults5 <- nbResults[!names(nbResults) %in% c("infectionDate.2")]
  expect_error(estimateRWrapper(nbResults5, indIDVar = "infectionDate"),
               "infectionDate.2 is not in the data frame.")
  
  #Testing timeFrame error
  expect_error(estimateRWrapper(nbResults, timeFrame = "garbage"),
               paste0("timeFrame must be one of: ",
                      paste0( c("days", "months", "weeks", "years"), collapse = ", ")))
})



test_that("Message printed in no range for average", {
  
  expect_message(estimateRWrapper(nbResults, rangeForAvg = NULL),
                 "Please choose the stable portion of the outbreak to calculate the average Rt")

  expect_message(estimateRWrapper(nbResults, rangeForAvg = NULL, bootSamples = 2),
                 "Please choose the stable portion of the outbreak to calculate the average Rt")
})


