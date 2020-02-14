
context("indToPair Function")
library(nbTransmission)

#Running on indData provided in package
pairU <- indToPair(indData = indData, indIDVar = "individualID", ordered = FALSE)
pairUD <- indToPair(indData = indData, indIDVar = "individualID",
                     dateVar = "infectionDate", units = "days", ordered = FALSE)
pairO <- indToPair(indData = indData, indIDVar = "individualID",
                   dateVar = "infectionDate", units = "hours", ordered = TRUE)

test_that("indToPair returns a dataframe of the correct length",{
  
  expect_true(is.data.frame(pairU))
  expect_true(is.data.frame(pairO))
  
  expect_true(nrow(pairU) == choose(nrow(indData), 2) * 2)
  expect_true(nrow(pairO) == choose(nrow(indData), 2))
})


test_that("indToPair with returns the right column names based on dateVar",{
  
  expect_true("infectionDate.Diff" %in% names(pairUD))
  expect_true("infectionDate.Diff" %in% names(pairO))
  
  expect_false("infectionDate.Diff" %in% names(pairU))
  
})


test_that("Descriptive error messages returned",{
  
  expect_error(indToPair(indData, indIDVar = "garbage"),
               "garbage is not in the data frame.")
  
  expect_error(indToPair(indData, indIDVar = "individualID", dateVar = "garbage"),
               "garbage is not in the data frame.")
  
  #Changing dates to character variables
  indData$infectionDatec <- as.character(indData$infectionDate)
  expect_error(indToPair(indData, indIDVar = "individualID", dateVar = "infectionDatec"),
               "infectionDatec must be either a date or a date-time (POSIXt) object.",
               fixed = TRUE)
  
  #Testing error about having ordered = TRUE with no date
  expect_error(indToPair(indData, indIDVar = "individualID", ordered = TRUE),
               "If ordered = TRUE, then dateVar must be provided")
  
  #Testing that missing units gives an error
  expect_error(indToPair(indData, indIDVar = "individualID", dateVar = "infectionDate"),
               "Please provide units for the time difference")
  
  #Providing an invalid clustering method
  expect_error(indToPair(indData, indIDVar = "individualID", dateVar = "infectionDate",
                         units = "garbage"),
               "units must be one of: mins, hours, days, weeks")
})


