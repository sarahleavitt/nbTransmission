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

pairData <- (covar[[1]]
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
pairData2 <- pairData %>% mutate(snpDist = ifelse(individualID.1 %in% trainingID |
                                          individualID.2 %in% trainingID, snpDist, NA))


#Adding these datasets to the package
use_data(indData, pairData2, overwrite = TRUE)

#Applying my function to get pairData
pairDataD <- indToPair(indData = indData, indIDVar = "individualID", separator = "_",
                       dateVar = "infectionDate", ordered = FALSE)



####################### Testing the Simulated Outbreak #######################


## Use the pairData dataset which represents a TB-like outbreak
# First create a dataset of ordered pairs
orderedPair <- pairData2[pairData2$infectionDiffY > 0, ]

## Create a variable called snpClose that will define probable links
# (<3 SNPs) and nonlinks (>12 SNPs) all pairs with between 2-12 SNPs
# will not be used to train.
orderedPair$snpClose <- ifelse(orderedPair$snpDist < 3, TRUE,
                               ifelse(orderedPair$snpDist > 12, FALSE, NA))
table(orderedPair$snpClose, useNA = "ifany")

## Running the algorithm
covariates = c("Z1", "Z2", "Z3", "Z4", "timeCat")
resGen <- calcProbabilities(orderedPair = orderedPair,
                            indIDVar = "individualID",
                            pairIDVar = "pairID",
                            goldStdVar = "snpClose",
                            covariates = covariates,
                            label = "SNPs", l = 1,
                            n = 10, m = 1, nReps = 10)

## Merging the probabilities back with the pair-level data
allProbs <- merge(resGen[[1]], orderedPair, by = "pairID", all = TRUE)

#Evaluating simulation
simEvaluate(allProbs)



## Calculating reproductive number ##

rInitial <- calcR(allProbs, dateVar = "infectionDate", indIDVar = "individualID",
                  pVar = "pScaled", timeFrame = "months")
rt <- rInitial[[2]]

#Cutting the outbreak
totalTime <- max(rt$timeRank, na.rm = TRUE) - min(rt$timeRank, na.rm = TRUE)
monthCut1 <- ceiling(0.1 * totalTime)
monthCut2 <- ceiling(0.7 * totalTime)

ggplot(data = rt, aes(x = timeRank, y = Rt)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(aes(xintercept = monthCut1), linetype = 2, size = 0.7) +
  geom_vline(aes(xintercept = monthCut2), linetype = 2, size = 0.7)

rFinal <- calcR(allProbs, dateVar = "infectionDate",
                     indIDVar = "individualID", pVar = "pScaled",
                     timeFrame = "months", rangeForAvg = c(monthCut1, monthCut2),
                     bootSamples = 100, alpha = 0.05)

rFinal[[3]]

ggplot(data = rFinal[[2]], aes(x = timeRank, y = Rt)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = ciLower, ymax = ciUpper), width = 0.1) +
  scale_y_continuous(name = "Monthly Effective Reproductive Number") + 
  scale_x_continuous(name = "Infection Year") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(aes(xintercept = monthCut1), linetype = 2, size = 0.7, color = "blue") +
  geom_vline(aes(xintercept = monthCut2), linetype = 2, size = 0.7, color = "blue") +
  geom_hline(data = rFinal[[3]], aes(yintercept = RtAvg), size = 0.7) +
  geom_hline(data = rFinal[[3]], aes(yintercept = ciLower), linetype = 2, 
             size = 0.7, color = "black") +
  geom_hline(data = rFinal[[3]], aes(yintercept = ciUpper), linetype = 2, 
             size = 0.7, color = "black")


## Creating Plots ##

nodes <- (indData
          %>% select(individualID, infectionDate)
          %>% arrange(infectionDate)
)
colBreaks <- c(-0.01, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)

#Clustering the infectors
clustRes <- (allProbs
             %>% group_by(individualID.2)
             %>% group_modify(~ findClustersKD(.x, pVar = "pScaled",
                                               binWidth = 0.04, minGap = 0))
)
table(clustRes$cluster)
topClust <- clustRes %>% filter(cluster == 1)
sum(topClust$transmission == TRUE) / length(unique(topClust$individualID.2))

edges <- (clustRes
          %>% select(individualID.1, individualID.2, transmission,
                     pScaled, snpDist, snpClose, cluster)
          #Arrange so high probability edges are drawn first
          %>% arrange(pScaled)
)
net <- graph_from_data_frame(d = edges, vertices = nodes, directed = T)
E(net)$pGroup <- cut(E(net)$pScaled, breaks = colBreaks, labels = 1:9)
#First get adjacency version of the network using pScaled
net.adj <- get.adjacency(net, attr = "pScaled", sparse = FALSE)


#Heatmap
#Marking cluster 1 with a *
net.trans <- get.adjacency(net, attr = "cluster", sparse = FALSE)
net.trans <- ifelse(net.trans == 1, "*", "")

par(mar = c(0, 0, 1, 0))
pheatmap(t(net.adj), cluster_rows = FALSE, cluster_cols = FALSE,
         border_color = NA, show_rownames = FALSE, show_colnames = FALSE,
         col=brewer.pal(9,"Blues"), breaks = colBreaks,
         display_numbers = t(net.trans), number_color = "white",
         fontsize_number = 7)

#Network of top cluster
par(mar = c(0, 0, 0.2, 0))
net_top <- delete.edges(net, E(net)[cluster == 2])
plot(net_top, vertex.size = 7, vertex.label.cex = 0.7,
     vertex.color = "gray",
     vertex.frame.color = "dark gray",
     edge.width = 2, edge.arrow.size = 0.4,
     edge.color = brewer.pal(9,"Blues")[E(net_top)$pGroup])


