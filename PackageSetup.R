
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

#Check for errors
devtools::check()


#### Running examples ####

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
edgeIDVar <- "edgeID"
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
resHam <- calcProbabilities(orderedPair = orderedHam, indIDVar = "individualID", edgeIDVar = "edgeID2",
                             goldStdVar = "SameGroup", covariates = covariates, label = "HamCont",
                             nbWeighting = FALSE, n = 10, m = 1, nReps = 5)

resHam2 <- full_join(orderedHam, resHam[[1]], by = "edgeID2")

rInitial <- calcR(resHam2, dateVar = "IsolationDate", indIDVar = "individualID",
                  pVar = "pScaled", timeFrame = "months")
rt <- rInitial[[2]]

#Cutting the outbreak
totalTime <- max(rt$timeRank) - min(rt$timeRank)
monthCut1 <- ceiling(0.1 * totalTime)
monthCut2 <- ceiling(0.9 * totalTime)

rFinal <- bootstrapR(resHam2, dateVar = "IsolationDate",
                     indIDVar = "individualID", pVar = "pScaled",
                     timeFrame = "months", rangeForAvg = c(monthCut1, monthCut2),
                     B = 1000, alpha = 0.05)

rFinal[[3]]


