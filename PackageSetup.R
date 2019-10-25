
setwd("~/Boston University/Dissertation/nbTransmission")
rm(list = ls())

library(devtools)
library(roxygen2)
library(qpdf)

# build()
# install()
# check()
# document()
# 
# #Add package dependence
# use_package()
# 
# #Add to ignore file
# use_build_ignore()
# 
# #Create a vingnette
# usethis::use_vignette("my-vignette")
# 
# #Create a directory of tests
# usethis::use_testthat()

#build, install, load
load_all()
document()

#Check for errors
devtools::check()


#### Starting tests ####
devtools::uses_testthat()
devtools::test()


#### Hamburg example ####

## Hamburg
hamInd <- readRDS("../Datasets/HamburgInd.rds")
hamPair <- readRDS("../Datasets/HamburgPair.rds")
hamPair <- hamPair %>% mutate(snpClose = ifelse(snpDist < 2, TRUE,
                                                ifelse(snpDist > 12, FALSE, NA)))

orderedHam <- (hamPair
               %>% mutate(IsolationDiff = as.numeric(difftime(IsolationDate.2,
                                                              IsolationDate.1, units = "days")))
               %>% filter(!is.na(IsolationDiff) & IsolationDiff > 0)
)

orderedPair <- orderedHam
dateVar <- "IsolationDate"
indIDVar <- "individualID"
pairIDVar <- "edgeID"
goldStdVar <- "snpClose"
pVar <- "pScaled"
label <- "Ham"
nbWeighting <- FALSE
n <- 10
m <- 1
nReps <- 1

covariates <- c("Study", "Nationality", "Sex", "Age", "SmearPos", "HIV",
                "SubstanceAbuse", "Residence", "Milieu", "TimeCat")

orderedHam <- orderedHam %>% rename(edgeID2 = edgeID)
resHam <- calcProbabilities(orderedPair = orderedHam, indIDVar = "individualID", pairIDVar = "edgeID2",
                             goldStdVar = "SameGroup", covariates = covariates, label = "HamCont",
                             n = 10, m = 1, nReps = 5)

resHam2 <- full_join(orderedHam, resHam[[1]], by = "edgeID2")

rInitial <- calcR(resHam2, dateVar = "IsolationDate", indIDVar = "individualID",
                  pVar = "pScaled", timeFrame = "months", bootSamples = 0)
rt <- rInitial[[2]]

#Cutting the outbreak
totalTime <- max(rt$timeRank) - min(rt$timeRank)
monthCut1 <- ceiling(0.1 * totalTime)
monthCut2 <- ceiling(0.9 * totalTime)

rFinal <- calcR(resHam2, dateVar = "IsolationDate",
                     indIDVar = "individualID", pVar = "pScaled",
                     timeFrame = "months", rangeForAvg = c(monthCut1, monthCut2),
                     bootSamples = 1000, alpha = 0.05)

rFinal[[3]]


#### Simulated Example ####

load_all()

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
resGen <- calcProbabilities(orderedPair = orderedPair,
                            indIDVar = "individualID",
                            pairIDVar = "pairID",
                            goldStdVar = "snpClose",
                            covariates = covariates,
                            label = "SNPs",
                            n = 10, m = 1, nReps = 1)

## Merging the probabilities back with the pair-level data
allProbs <- resGen[[1]] %>% full_join(orderedPair, by = "pairID")
summary(allProbs$pScaled)



rInitial <- calcR(allProbs, dateVar = "infectionDate", indIDVar = "individualID",
                  pVar = "pScaled", timeFrame = "months", bootSamples = 0)
rt <- rInitial[[2]]
names(rInitial)
names(rInitial$RiDf)
names(rInitial$RtDf)
names(rInitial$RtAvgDf)

#Cutting the outbreak
totalTime <- max(rt$timeRank) - min(rt$timeRank)
monthCut1 <- ceiling(0.1 * totalTime)
monthCut2 <- ceiling(0.7 * totalTime)

rFinal <- calcR(allProbs, dateVar = "infectionDate", indIDVar = "individualID",
                pVar = "pScaled", timeFrame = "months",
                rangeForAvg = c(monthCut1, monthCut2),
                bootSamples = 1000, alpha = 0.05)
rFinal[[3]]

names(rFinal)
names(rFinal$RiDf)
names(rFinal$RtDf)
names(rFinal$RtAvgDf)

