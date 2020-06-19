
context("visualizeResults Functions")
library(nbTransmission)


#### Testing Network and Heatmap Plots ####

nbNetworkWrapper <- function(nbResults,
                             dateVar = "infectionDate",
                             indIDVar = "individualID",
                             pVar = "pScaled",
                             clustMethod = "none",
                             cutoff = NA,
                             blackAndWhite = FALSE,
                             probBreaks = c(-0.01, 0.001, 0.005, 0.01,
                                            0.05, 0.1, 0.25, 0.5, 0.75, 1)){
  
  nbNetwork(nbResults, indIDVar = indIDVar, dateVar = dateVar,
          pVar = pVar, clustMethod = clustMethod, cutoff = cutoff,
          blackAndWhite = blackAndWhite, probBreaks = probBreaks)
}

nbHeatmapWrapper <- function(nbResults,
                             dateVar = "infectionDate",
                             indIDVar = "individualID",
                             pVar = "pScaled",
                             clustMethod = "none",
                             cutoff = NA,
                             blackAndWhite = FALSE,
                             probBreaks = c(-0.01, 0.001, 0.005, 0.01,
                                            0.05, 0.1, 0.25, 0.5, 0.75, 1)){
  
  nbHeatmap(nbResults, indIDVar = indIDVar, dateVar = dateVar,
            pVar = pVar, clustMethod = clustMethod, cutoff = cutoff,
            blackAndWhite = blackAndWhite, probBreaks = probBreaks)
}

#Internal function for both nbNetwork and nbHeatmap
createNetworkWrapper <- function(nbResults,
                             dateVar = "infectionDate",
                             indIDVar = "individualID",
                             pVar = "pScaled",
                             clustMethod = "none",
                             cutoff = NA,
                             probBreaks = c(-0.01, 0.001, 0.005, 0.01,
                                            0.05, 0.1, 0.25, 0.5, 0.75, 1)){
  
  createNetwork(nbResults, indIDVar = indIDVar, dateVar = dateVar,
            pVar = pVar, clustMethod = clustMethod, cutoff = cutoff,
            probBreaks = probBreaks)
}

#Shortening dataset for sake of speed
testData <- nbResults[1:500, ]


test_that("Plot functions return null objects with no errors",{
  
  net1 <- nbNetworkWrapper(testData)
  net2 <- nbNetworkWrapper(testData, clustMethod = "hc_absolute", cutoff = 0.05)
  net3 <- nbNetworkWrapper(testData, clustMethod = "hc_absolute", cutoff = 0.05,
                           blackAndWhite = TRUE)
  
  expect_true(class(net1) == "NULL")
  expect_true(class(net2) == "NULL")
  expect_true(class(net3) == "NULL")
  
  heat1 <- nbHeatmapWrapper(testData)
  heat2 <- nbHeatmapWrapper(testData, clustMethod = "hc_absolute", cutoff = 0.05)
  heat3 <- nbHeatmapWrapper(testData, clustMethod = "hc_absolute", cutoff = 0.05,
                            blackAndWhite = TRUE)
  
  expect_true(class(heat1) == "pheatmap")
  expect_true(class(heat2) == "pheatmap")
  expect_true(class(heat3) == "pheatmap")
  
})


test_that("Internal createNetwork function returns igraph object",{
  
  net1 <- createNetworkWrapper(testData)
  net2 <- createNetworkWrapper(testData, clustMethod = "hc_absolute", cutoff = 0.05)
  
  expect_true(class(net1) == "igraph")
  expect_true(class(net2) == "igraph")
  
})


#Only need to test for the internal function because both nbNetwork and nbHeatmap call it
test_that("Descriptive error messages returned from internal createNetwork function",{
  
  expect_error(createNetworkWrapper(testData, indIDVar = "garbage"),
               "garbage.1 is not in the data frame.")
  
  expect_error(createNetworkWrapper(testData, dateVar = "garbage"),
               "garbage.1 is not in the data frame.")
  
  expect_error(createNetworkWrapper(testData, pVar = "garbage"),
               "garbage is not in the data frame.")
  
  #Removing individualID columns
  testData2 <- testData[!names(testData) %in% c("individualID.1")]
  expect_error(createNetworkWrapper(testData2, indIDVar = "individualID"),
               "individualID.1 is not in the data frame.")
  
  testData3 <- testData[!names(testData) %in% c("individualID.2")]
  expect_error(createNetworkWrapper(testData3, indIDVar = "individualID"),
               "individualID.2 is not in the data frame.")
  
  #Removing the date columns
  testData4 <- testData[!names(testData) %in% c("infectionDate.1")]
  expect_error(createNetworkWrapper(testData4, dateVar = "infectionDate"),
               "infectionDate.1 is not in the data frame.")
  
  testData5 <- testData[!names(testData) %in% c("infectionDate.2")]
  expect_error(createNetworkWrapper(testData5, dateVar = "infectionDate"),
               "infectionDate.2 is not in the data frame.")
  
  #Changing dates to character variables
  testData$infectionDatec.1 <- as.character(testData$infectionDate.1)
  testData$infectionDatec.2 <- testData$infectionDate.2
  expect_error(createNetworkWrapper(testData, dateVar = "infectionDatec"),
               "infectionDatec.1 must be either a date or a date-time (POSIXt) object.",
               fixed = TRUE)
  
  testData$infectionDatec.1 <- testData$infectionDate.1
  testData$infectionDatec.2 <- as.character(testData$infectionDate.2)
  expect_error(createNetworkWrapper(testData, dateVar = "infectionDatec"),
               "infectionDatec.2 must be either a date or a date-time (POSIXt) object.",
               fixed = TRUE)
  
  testData5 <- testData[!names(testData) %in% c("infectionDate.2")]
  expect_error(createNetworkWrapper(testData5, indIDVar = "infectionDate"),
               "infectionDate.2 is not in the data frame.")
  
  
  #Testing that misclustng clustMethod gets set to none and gives warning
  expect_warning(createNetwork(testData, dateVar = "infectionDate", indIDVar = "individualID",
                           pVar = "pScaled"),
                 "No clustMethod was provided so it was set to 'none'")
  
  #Providing an invalid clustering method
  expect_error(createNetworkWrapper(testData, clustMethod = "garbage"),
               "clustMethod must be one of: none, n, kd, hc_absolute, hc_relative")
  
  #Providing a clust method with no cutoff
  expect_error(createNetworkWrapper(testData, clustMethod = "hc_absolute"),
               "Please provide one or more cutoff values")

})


#Only need to test for the internal function because both nbNetwork and nbHeatmap call it
test_that("Descriptive error messages returned for probBreaks",{
  
  expect_error(createNetworkWrapper(nbResults, probBreaks = c(-0.01, 5, 1)),
               "All values of probBreaks should be less than 1")
  
  expect_error(createNetworkWrapper(nbResults, probBreaks = c(-0.01, 1)),
               "Please make sure probBreaks has between 3 and 10 elements")
  
  expect_error(createNetworkWrapper(nbResults, probBreaks = c(-0.01, seq(0.01, 0.99, 0.01), 1)),
               "Please make sure probBreaks has between 3 and 10 elements")
  
  expect_message(createNetworkWrapper(nbResults, probBreaks = c(0.01, 0.05, 1)),
                 "First element of probBreaks is not negative so -0.01 was added to the beginning")
  
  expect_message(createNetworkWrapper(nbResults, probBreaks = c(-0.01, 0.05)),
                 "Last element of probBreaks is not 1 so 1 was added to the end")
  
  expect_message(createNetworkWrapper(nbResults, probBreaks = c(0.01, 0.05)),
                 "First element of probBreaks is not negative so -0.01 was added to the beginning")
  
  expect_message(createNetworkWrapper(nbResults, probBreaks = c(0.01, 0.05)),
                 "Last element of probBreaks is not 1 so 1 was added to the end")
  
})



#### Testing Rt Plot ####

rFinal <- estimateR(nbResults, dateVar = "infectionDate",
                    indIDVar = "individualID", pVar = "pScaled",
                    timeFrame = "months", rangeForAvg = c(10, 150),
                    bootSamples = 2, alpha = 0.05)

rFinal2 <- estimateR(nbResults, dateVar = "infectionDate",
                    indIDVar = "individualID", pVar = "pScaled",
                    timeFrame = "days", alpha = 0.05)


test_that("Plot functions return null objects with no errors",{
  
  rt1 <- plotRt(rFinal, includeRtAvg = TRUE, includeRtCI = TRUE, includeRtAvgCI = TRUE)
  rt2 <- plotRt(rFinal, includeRtAvg = TRUE, includeRtCI = TRUE, includeRtAvgCI = FALSE)
  rt3 <- plotRt(rFinal, includeRtAvg = FALSE, includeRtCI = TRUE, includeRtAvgCI = FALSE)
  rt4 <- plotRt(rFinal, includeRtAvg = TRUE, includeRtCI = FALSE, includeRtAvgCI = TRUE)
  rt5 <- plotRt(rFinal, includeRtAvg = TRUE, includeRtCI = FALSE, includeRtAvgCI = FALSE)
  rt6 <- plotRt(rFinal, includeRtAvg = FALSE, includeRtCI = FALSE, includeRtAvgCI = FALSE)
  rt7 <- plotRt(rFinal2, includeRtAvg = TRUE)
  
  expect_true("ggplot" %in% class(rt1))
  expect_true("ggplot" %in% class(rt2))
  expect_true("ggplot" %in% class(rt3))
  expect_true("ggplot" %in% class(rt4))
  expect_true("ggplot" %in% class(rt5))
  expect_true("ggplot" %in% class(rt6))
  expect_true("ggplot" %in% class(rt7))
  
})


test_that("Descriptive error messages returned for plotRt", {
  
  expect_error(plotRt("garbage"),
               "The rData argument should be the list output from the function estimateR")
  
  expect_error(plotRt(list("garbage")),
               "The rData argument should be the list output from the function estimateR")
  
  expect_error(plotRt(rFinal2, includeRtAvgCI = TRUE),
               "Please provide a rData list that has confidence intervals")
  
  expect_error(plotRt(rFinal2, includeRtCI = TRUE),
               "Please provide a rData list that has confidence intervals")

  expect_error(plotRt(rFinal2, includeRtCI = TRUE, includeRtAvgCI = TRUE),
               "Please provide a rData list that has confidence intervals")
})


