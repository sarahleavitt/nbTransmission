
setwd("~/Boston University/Dissertation/nbTransmission")
rm(list = ls())

library(devtools)
library(roxygen2)

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



#### Running examples ####

## Hamburg
hamInd <- readRDS("../Datasets/HamburgInd.rds")
hamPair <- readRDS("../Datasets/HamburgPair.rds")
hamPair <- hamPair %>% mutate(snpClose = ifelse(snpDist < 2, TRUE,
                                         ifelse(snpDist > 12, FALSE, NA)))

orderedHam <- (hamPair
               %>% mutate(IsolationDiff = as.numeric(difftime(IsolationDate.2,
                                                              IsolationDate.1, units = "days")))
               %>% filter(!is.na(IsolationDiff) & observationDiff > 0)
)

# orderedPair <- orderedHam
# dateVar <- "IsolationDate"
# indIDVar <- "individualID"
# edgeIDVar <- "edgeID"
# goldStdVar <- "snpClose"
# pVar <- "pScaled"
# label <- "Ham"
# nbWeighting <- FALSE
# n <- 10
# m <- 1
# nReps <- 1

load_all()

set.seed(103020)
covariates <- c("Study", "Nationality", "Sex", "Age", "SmearPos", "HIV",
                "SubstanceAbuse", "Residence", "Milieu", "TimeCat")

resHam1 <- calcProbabilities(orderedPair = orderedHam, indIDVar = "individualID", edgeIDVar = "edgeID",
                              goldStdVar = "snpClose", covariates = covariates, label = "Ham",
                              nbWeighting = FALSE, n = 10, m = 1, nReps = 50)

resHam2 <- calcProbabilities(orderedPair = orderedHam, indIDVar = "individualID", edgeIDVar = "edgeID",
                              goldStdVar = "SameGroup", covariates = covariates, label = "Ham",
                              nbWeighting = FALSE, n = 10, m = 1, nReps = 50)

resHamCov1 <- resHam1[[1]] %>% full_join(orderedHam, by = "edgeID")
resHamCov2 <- resHam2[[1]] %>% full_join(orderedHam, by = "edgeID")
hamRes <- bind_rows(resHamCov1, resHamCov2)

ri1 <- calcRi(hamRes1, dateVar = "IsolationDate", indIDVar = "individualID", pVar = "pScaled")
ri2 <- calcRi(hamRes2, dateVar = "IsolationDate", indIDVar = "individualID", pVar = "pScaled")

riSum <- (hamRes
          %>% group_by(label)
          %>% do(calcRi(., dateVar = "IsolationDate", indIDVar = "individualID", pVar = "pScaled"))
)

#Calculating the average monthly Rt by run
repNumMH <- (riSum
             %>% group_by(label, monthR)
             %>% summarize(Rt = mean(Ri),
                           month = first(month),
                           year = first(year))
             %>% ungroup()
)

#Cutting the outbreak
totalTime <- max(repNumMH$monthR) - min(repNumMH$monthR)
monthCut1 <- ceiling(0.1 * totalTime)
monthCut2 <- ceiling(0.9 * totalTime)

repNumCutH <- repNumMH  %>% filter(monthR > monthCut1 & monthR < monthCut2)

#Calculating the average R0 for each method
monthR0H <- (repNumCutH
             %>% group_by(label)
             %>% summarize(R0 = mean(Rt, na.rm = TRUE),
                           R0_median = median(Rt, na.rm = TRUE),
                           sdR0 = sd(Rt, na.rm = TRUE))
)
monthR0H
tapply(hamRes$pAvg, hamRes$label, summary)


## Simulation
sampleSize <- 200
neg <- 0.25
pi <- 1
off.r <- 1.2
off.p <- 0.5
w.shape <- 1.05
w.scale <- 1 / (0.0014 * 365)
ws.shape <- w.shape
ws.scale <- w.scale
shift <- 0.25
multOutbreaks <- TRUE
rootseq <- NULL
length <- 300
rate <- 0.5 / length

source("../dissertation_code/SimOutbreak.R")
source("../dissertation_code/SimCovariates.R")
source("../dissertation_code/TransPhylo/SimulateOutbreakS.R")

#Simulate outbreak  
set.seed(1001)
obk <- simOutbreak(neg = neg, pi = pi, off.r = off.r, off.p = off.p,
                   w.scale = w.scale, w.shape = w.shape, shift = shift,
                   ws.scale = ws.scale, ws.shape = ws.shape,
                   sampleSize = sampleSize,
                   multOutbreaks = multOutbreaks,
                   length = length, rate = rate)
initialInd <- obk[[1]]
initialPair <- obk[[2]]

#Simulating covariates
covar <- simCovariates(initialInd, initialPair, scheme = "info")
indData <- covar[[2]]
pairData <- covar[[1]]

#Subseting to the pairs with the potential infector observed before the infectee
#Restricting to sampled cases
orderedPair <- pairData %>% filter(!is.na(observationDiff) & observationDiff > 0)


dateVar <- "infectionDate"
indIDVar <- "individualID"
edgeIDVar <- "edgeID"
goldStdVar <- "transmission"
pVar <- "pScaled"
covariates <- c("Y1", "Y2", "Y3", "Y4", "timeCat")

results2 <- calcProbabilities(orderedPair, indIDVar, edgeIDVar, goldStdVar, 
                              covariates, label = NULL, n = 10, m = 1, nReps = 1)

df <- results2[[1]] %>% full_join(orderedPair, by = "edgeID")

ri2 <- calcRi(df, dateVar, indIDVar, pVar)


