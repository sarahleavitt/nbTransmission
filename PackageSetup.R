
setwd("~/Boston University/Dissertation/nbTransmission")

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

orderedPair <- orderedHam
dateVar <- "IsolationDate"
indIDVar <- "individualID"
edgeIDVar <- "edgeID"
goldStdVar <- "snpClose"
pVar <- "pScaled"
covariates <- c("Study", "Nationality", "Sex", "Age", "SmearPos", "HIV",
               "SubstanceAbuse", "Residence", "Milieu", "TimeCat")
label <- "Ham"
nbWeighting <- FALSE
n <- 10
m <- 1
nReps <- 1

load_all()
results1 <- calcProbabilities(orderedHam, indIDVar, edgeIDVar, goldStdVar,
                              covariates, label = "Ham", nbWeighting = FALSE,
                              n = 10, m = 1, nReps = 1)

hamRes <- results1[[1]] %>% full_join(orderedHam, by = "edgeID")

ri1 <- calcRi(hamRes, dateVar, indIDVar, pVar)


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


