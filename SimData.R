
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
library(glmnet)
library(caret)
library(gtools)

library(igraph)
library(RColorBrewer)
library(pheatmap)

################ Creating a Simulated Outbreak ###############

source("../dissertation_code/SimOutbreak.R")
source("../dissertation_code/SimulateOutbreakS.R")
source("../dissertation_code/SimCovariates.R")
source("../dissertation_code/SimEvaluate.R")
source("../dissertation_code/CalcSI.R")
load_all()

neg <- 0.25
pi <- 1
off.r <- 1.2
off.p <- 0.5
w.shape <- 1.2
w.scale <- 2
ws.shape <- w.shape
ws.scale <- w.scale
shift <- 0
multOutbreaks <- FALSE 
rootseq <- NULL
length <- 300
rate <- 0.5 / length
time <- 20
sampleSize <- 100

#Simulate outbreak  
set.seed(1004)
obk <- simOutbreak(neg = neg, pi = pi, off.r = off.r, off.p = off.p,
                   w.scale = w.scale, w.shape = w.shape, shift = shift,
                   ws.scale = ws.scale, ws.shape = ws.shape,
                   sampleSize = sampleSize, time = time,
                   multOutbreaks = multOutbreaks,
                   length = length, rate = rate)
print(paste0("Simulated outbreak, n = ", nrow(obk[[1]])))

#Simulating covariates
covar <- simCovariates(obk[[1]], obk[[2]], scheme = "info")
print("Simulated covariates")

indData <- (covar[[2]]
            %>% filter(!is.na(sampleDate))
            %>% select(individualID, infector, infectionDate, sampleDate,
                                X1, X2, X3, X4)
            %>% mutate(individualID = as.numeric(gsub("^10", "", individualID)),
                       infector = as.numeric(gsub("^10", "", infector)))
)
nrow(indData)
table(indData$individualID)

pairData1 <- (covar[[1]]
             %>% filter(!is.na(sampleDate.1) & !is.na(sampleDate.2))
             %>% mutate(infectionDiffY = infectionDiff / 365,
                        individualID.1 = as.numeric(gsub("^10", "", individualID.1)),
                        individualID.2 = as.numeric(gsub("^10", "", individualID.2)),
                        pairID = gsub("_10", "_", gsub("^10", "", edgeID)))
             %>% select(pairID, individualID.1, individualID.2, transmission, snpDist,
                        infectionDate.1, infectionDate.2, sampleDate.1, sampleDate.2,
                        sampleDiff, infectionDiff, infectionDiffY, timeCat,
                        Z1 = Y1, Z2 = Y2, Z3 = Y3, Z4 = Y4)
)

#Removing snpDist from a random set of individuals to better replicate real data
pTraining <- 0.6
trainingID <- (indData
               %>% sample_frac(pTraining)
               %>% pull(individualID)
)
pairData <- pairData1 %>% mutate(snpDist = ifelse(individualID.1 %in% trainingID |
                                          individualID.2 %in% trainingID, snpDist, NA))

#Adding these datasets to the package
#use_data(indData, pairData, overwrite = TRUE)


#Applying my function to get pairData
pairDataD <- indToPair(indData = indData, indIDVar = "individualID", separator = "_",
                       dateVar = "infectionDate", ordered = FALSE)



####################### Testing the Simulated Outbreak #######################


## Use the pairData dataset which represents a TB-like outbreak
# First create a dataset of ordered pairs
orderedPair <- pairData[pairData$infectionDiffY > 0, ]

## Create a variable called snpClose that will define probable links
# (<3 SNPs) and nonlinks (>12 SNPs) all pairs with between 2-12 SNPs
# will not be used to train.
orderedPair$snpClose <- ifelse(orderedPair$snpDist < 3, TRUE,
                               ifelse(orderedPair$snpDist > 10, FALSE, NA))
table(orderedPair$snpClose, useNA = "ifany")

## Running the algorithm
set.seed(0)
covariates = c("Z1", "Z2", "Z3", "Z4", "timeCat")
resGen <- nbProbabilities(orderedPair = orderedPair,
                            indIDVar = "individualID",
                            pairIDVar = "pairID",
                            goldStdVar = "snpClose",
                            covariates = covariates,
                            label = "SNPs", l = 1,
                            n = 10, m = 1, nReps = 10)

## Merging the probabilities back with the pair-level data
nbResults <- merge(resGen[[1]], orderedPair, by = "pairID", all = TRUE)

#Evaluating simulation
simEvaluate(nbResults)

#use_data(nbResults, overwrite = TRUE)



## Calculating reproductive number ##

rInitial <- estimateR(nbResults, dateVar = "infectionDate", indIDVar = "individualID",
                  pVar = "pScaled", timeFrame = "months")
rt <- rInitial$RtDf

plotRt(rInitial, includeRtAvg = FALSE, includeRtCI = FALSE, includeRtAvgCI = FALSE)

#Cutting the outbreak
totalTime <- max(rt$timeRank, na.rm = TRUE) - min(rt$timeRank, na.rm = TRUE)
cut1 <- ceiling(0.15 * totalTime)
cut2 <- ceiling(0.85 * totalTime)

cut1 <- 25
cut2 <- 125

ggplot(data = rt, aes(x = timeRank, y = Rt)) +
  geom_point() +
  geom_line() +
  geom_hline(data = rInitial$RtAvgDf, aes(yintercept = RtAvg), size = 0.7) +
  geom_vline(aes(xintercept = cut1), linetype = 2, size = 0.7) +
  geom_vline(aes(xintercept = cut2), linetype = 2, size = 0.7)

rFinal <- estimateR(nbResults, dateVar = "infectionDate",
                     indIDVar = "individualID", pVar = "pScaled",
                     timeFrame = "months", rangeForAvg = c(cut1, cut2),
                     bootSamples = 100, alpha = 0.05)

rFinal$RtAvgDf

#Plot of reproductive number over time
ggplot(data = rFinal[[2]], aes(x = timeRank, y = Rt)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = ciLower, ymax = ciUpper), width = 0.1) +
  scale_y_continuous(name = "Monthly Effective Reproductive Number") + 
  scale_x_continuous(name = "Infection Month") +
  geom_vline(aes(xintercept = cut1), linetype = 2, size = 0.7, color = "blue") +
  geom_vline(aes(xintercept = cut2), linetype = 2, size = 0.7, color = "blue") +
  geom_hline(data = rFinal[[3]], aes(yintercept = RtAvg), size = 0.7) +
  geom_hline(data = rFinal[[3]], aes(yintercept = ciLower), linetype = 2, 
             size = 0.7, color = "black") +
  geom_hline(data = rFinal[[3]], aes(yintercept = ciUpper), linetype = 2, 
             size = 0.7, color = "black")


plotRt(rInitial, includeRtAvg = FALSE, includeRtCI = FALSE, includeRtAvgCI = FALSE)
plotRt(rFinal, includeRtAvg = TRUE, includeRtCI = TRUE, includeRtAvgCI = TRUE)




## Clustering Probabilities ##

## Clustering using top n
# Top cluster includes infectors with highest 3 probabilities
clust1 <- clusterInfectors(nbResults, indIDVar = "individualID", pVar = "pScaled",
                           clustMethod = "n", cutoff = 3)
table(clust1$cluster)

## Clustering using hierarchical clustering

# Cluster all infectees, do not force gap to be certain size
clust2 <- clusterInfectors(nbResults, indIDVar = "individualID", pVar = "pScaled",
                           clustMethod = "hc_absolute", cutoff = 0)
table(clust2$cluster)

# Absolute difference: gap between top and bottom clusters is more than 0.05
clust3 <- clusterInfectors(nbResults, indIDVar = "individualID", pVar = "pScaled",
                           clustMethod = "hc_absolute", cutoff = 0.05)
table(clust3$cluster)

# Relative difference: gap between top and bottom clusters is more than double any other gap
clust4 <- clusterInfectors(nbResults, indIDVar = "individualID", pVar = "pScaled",
                           clustMethod = "hc_relative", cutoff = 2)
table(clust4$cluster)

## Clustering using kernel density estimation
# Using a small binwidth of 0.01
clust5 <- clusterInfectors(nbResults, indIDVar = "individualID", pVar = "pScaled",
                           clustMethod = "kd", cutoff = 0.01)
table(clust5$cluster)


## Estimating the Serial Interval ##

#Rerunning without time
resGenNoT <- nbProbabilities(orderedPair = orderedPair,
                            indIDVar = "individualID",
                            pairIDVar = "pairID",
                             goldStdVar = "snpClose",
                             covariates = c("Z1", "Z2", "Z3", "Z4"),
                             label = "SNPs", l = 1,
                             n = 10, m = 1, nReps = 1)

nbResultsNoT <- merge(resGen[[1]], orderedPair, by = "pairID", all = TRUE)

dateVar <- "infectionDiffY"
#Finding the true values
truePairs <- nbResults[nbResults$transmission == TRUE, ]
obsMean <- mean(truePairs[, dateVar])
obsMedian <- median(truePairs[, dateVar])
obsSD <- sd(truePairs[, dateVar])

#Using wrapper
estimateSI(nbResults, indIDVar = "individualID", timeDiffVar = "infectionDiffY",
           pVar = "pScaled", clustMethod = "none", initialPars = c(2, 2))

siPars <- estimateSI(nbResults, indIDVar = "individualID", timeDiffVar = "infectionDiffY",
           pVar = "pScaled", clustMethod = "hc_absolute", cutoff = 0.05,
           initialPars = c(2, 2), bootSamples = 0)

siParsCI <- estimateSI(nbResults, indIDVar = "individualID", timeDiffVar = "infectionDiffY",
           pVar = "pScaled", clustMethod = "hc_absolute", cutoff = 0.05,
           initialPars = c(2, 2), bootSamples = 5)

trainPairs <- nbResults[nbResults$snpClose == TRUE, ]
siParsGS <- estimateSI(trainPairs, indIDVar = "individualID", timeDiffVar = "infectionDiffY",
                       pVar = "pScaled", clustMethod = "none", initialPars = c(2, 2))

#Plotting the results
ggplot(data = nbResults, aes(x = infectionDiffY)) +
  geom_histogram(data = truePairs, aes(y = ..density..), bins = 40) +
  scale_x_continuous(name = "Serial Interval (years)", limits = c(0, 20)) +
  geom_line(aes(y = dgammaS(infectionDiffY, shape = siPars$shape, scale = siPars$scale,
                            shift = 0), color = "Top Cluster"), size = 1.2) +
  geom_line(aes(y = dgammaS(infectionDiffY, shape = siParsCI$shapeCILB,
                            scale = siParsCI$scaleCIUB,
                            shift = 0), color = "CI Lower"), size = 1.2) +
  geom_line(aes(y = dgammaS(infectionDiffY, shape = siParsCI$shapeCIUB,
                            scale = siParsCI$scaleCIUB,
                            shift = 0), color = "CI Upper"), size = 1.2)
  

#Using a shift
estimateSIPars(nbResults, indIDVar = "individualID", timeDiffVar = "infectionDiffY",
           pVar = "pScaled", clustMethod = "hc_absolute", cutoff = 0.05,
           initialPars = c(2, 2), shift = 0.25)


#Using clusterInfectors and performPEM
performPEM(nbResults, indIDVar = "individualID", timeDiffVar = "infectionDiffY",
           pVar = "pScaled", initialPars = c(2, 2), shift = 0, plot = TRUE)

allClust <- clusterInfectors(nbResults, indIDVar = "individualID", pVar = "pScaled",
                             clustMethod = "hc_absolute", cutoff = 0.05)

performPEM(allClust[allClust$cluster == 1, ], indIDVar = "individualID",
           timeDiffVar = "infectionDiffY", pVar = "pScaled",
           initialPars = c(2, 2), shift = 0, plot = TRUE)



## Creating Plots ##


nbHeatmap(nbResults, indIDVar = "individualID", dateVar = "infectionDate",
          pVar = "pScaled", clustMethod = "none")

nbHeatmap(nbResults, indIDVar = "individualID", dateVar = "infectionDate",
          pVar = "pScaled", clustMethod = "hc_absolute", cutoff = 0.05,
          blackAndWhite = TRUE, probBreaks = c(-0.01, 0.01, 0.1, 0.25, 0.5, 1))

## Network of all pairs in color with the default probability breaks
nbNetwork(nbResults, indIDVar = "individualID", dateVar = "infectionDate",
pVar = "pScaled", clustMethod = "none")

## Adding stars for the top cluster, in black and white, changing the probability breaks
nbNetwork(nbResults, indIDVar = "individualID", dateVar = "infectionDate",
          pVar = "pScaled", clustMethod = "hc_absolute", cutoff = 0.05,
          blackAndWhite = TRUE, probBreaks = c(-0.01, 0.01, 0.1, 0.25, 0.5, 1))


## Plotting probabilities ##

ggplot(data = nbResults %>% filter(individualID.2 < 50),
       aes(x = pRank, y = pScaled, color = transmission)) +
  geom_point() +
  facet_wrap(~individualID.2) +
  theme(legend.position = "none")

ggplot(data = nbResults %>% filter(individualID.2 > 50),
       aes(x = pRank, y = pScaled, color = transmission)) +
  geom_point() +
  facet_wrap(~individualID.2) +
  theme(legend.position = "none")


