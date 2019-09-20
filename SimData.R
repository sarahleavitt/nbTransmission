## Created Simulated Dataset ##

setwd("~/Boston University/Dissertation/nbTransmission")
rm(list = ls())

library(devtools)
library(roxygen2)
library(qpdf)
library(TransPhylo)
library(phangorn)
library(reshape2)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)

################ Simulate Outbreak ###############

source("../dissertation_code/SimOutbreak.R")
source("../dissertation_code/SimulateOutbreakS.R")
source("../dissertation_code/SimCovariates.R")

neg <- 0.25
pi <- 1
off.r <- 1.2
off.p <- 0.5
w.shape <- 1.2
w.scale <- 2
ws.shape <- w.shape
ws.scale <- w.scale
shift <- 0.25
multOutbreaks <- TRUE 
rootseq <- NULL
length <- 300
rate <- 0.5 / length

sampleSize <- 100

#Simulate outbreak  
set.seed(10010)
obk <- simOutbreak(neg = neg, pi = pi, off.r = off.r, off.p = off.p,
                   w.scale = w.scale, w.shape = w.shape, shift = shift,
                   ws.scale = ws.scale, ws.shape = ws.shape,
                   sampleSize = sampleSize, multOutbreaks = TRUE,
                   length = length, rate = rate)
print(paste0("Simulated outbreak, n = ", nrow(indData)))

#Simulating covariates
covar <- simCovariates(obk[[1]], obk[[2]], scheme = "info")
pairData <- covar[[1]]
indData <- covar[[2]]
print("Simulated covariates")

indData <- indData %>% select(individualID, outbreakID, infector, infectionDate, sampleDate,
                                X1, X2, X3, X4)
pairData <- (pairData
             %>% select(-observationDate.1, -observationDate.2, -observationDiff,
                                  -Y1.D, -Y2.D, -Y3.D, -Y4.D, -Y1.C, -Y2.C, -Y3.C, -Y4.C,
                                  Z1 = Y1, Z2 = Y2, Z3 = Y3, Z4 = Y4, -complete)
             %>% mutate(infectionDiffY = infectionDiff / 365)
)

#Adding these datasets to the package
use_data(indData, pairData, overwrite = TRUE)

pTraining <- 1
#Finding all pairs that can be included in the training dataset (pTraining)
trainingID <- (indData
               %>% filter(!is.na(sampleDate))
               %>% sample_frac(pTraining)
               %>% pull(individualID)
)

orderedPair <- (pairData
                %>% filter(infectionDiffY > 0)
                %>% mutate(snpClose = ifelse(snpDist < 2, TRUE,
                                      ifelse(snpDist > 12, FALSE, NA)),
                           trainPair = ifelse(individualID.1 %in% trainingID &
                                                individualID.2 %in% trainingID &
                                                !is.na(snpClose),
                                              TRUE, FALSE),
                           snpCloseGS = ifelse(trainPair == TRUE, snpClose, NA))
)

orderedPair <- pairData[pairData$infectionDiffY > 0, ]
orderedPair$snpClose <- ifelse(orderedPair$snpDist < 3, TRUE,
                               ifelse(orderedPair$snpDist > 12, FALSE, NA))
table(orderedPair$snpClose)
table(orderedPair$transmission, orderedPair$snpClose)
prop.table(table(orderedPair$transmission, orderedPair$snpClose, useNA = "ifany"), 1)

load_all()
covariates = c("Z1", "Z2", "Z3", "Z4", "timeCat")
resGen <- calcProbabilities(orderedPair = orderedPair, indIDVar = "individualID",
                            edgeIDVar = "edgeID", goldStdVar = "snpClose", covariates = covariates,
                            label = "SNPs", nbWeighting = FALSE, n = 10, m = 1, nReps = 10)
allProbs <- resGen[[1]] %>% full_join(orderedPair, by = "edgeID")


rInitial <- calcR(allProbs, dateVar = "infectionDate", indIDVar = "individualID",
                  pVar = "pScaled", timeFrame = "months")
rt <- rInitial[[2]]


#Cutting the outbreak
totalTime <- max(rt$timeRank) - min(rt$timeRank)
monthCut1 <- ceiling(0.1 * totalTime)
monthCut2 <- ceiling(0.7 * totalTime)

rFinal <- bootstrapR(allProbs, dateVar = "infectionDate",
                     indIDVar = "individualID", pVar = "pScaled",
                     timeFrame = "months", rangeForAvg = c(monthCut1, monthCut2),
                     B = 10, alpha = 0.05)

rFinal[[3]]


