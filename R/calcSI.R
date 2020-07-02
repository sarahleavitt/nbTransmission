
#' Estimates the generation/serial interval distribution
#'
#' The function \code{estimateSI} uses the relative transmission probabilities to estimate
#' the generation/serial interval distribution assuming a gamma distribution.
#' It uses the PEM algorithm developed by Hens et al. 2012 extending their method to include
#' restricting analysis to the top cluster of possible infectors.
#' 
#' The PEM algorithm uses the prior probability that each pair is connected by direct
#' transmission to estimate the generation/serial interval using estimation maximization.
#' This code will provide an estimate of the generation interval if \code{timeDiffVar} represents the
#' difference in infection dates and the serial interval if it represents the difference in symptom onset dates.
#' The current generation/serial interval distribution parameters are used to update the probabilities
#' which are then used to update the generation/serial interval distribution parameters. The process
#' continues until the parameters converge (indicated by a change of less than \code{epsilon})
#' between iterations. \emph{Note: time difference between pairs should not be used to 
#' estimate the probabilities}
#' 
#' This function acts as a wrapper around \code{\link{performPEM}} which integrates
#' estimation of the generation/serial interval distribution with clustering the infectors and calculates
#' derived parameters (mean, median, sd) of the distribution. Generally, this function
#' should be called instead of \code{\link{performPEM}} directly.
#' 
#' All pairs of cases can be used in the estimation process by setting
#' \code{clustMethod = "none"}. However, if the probabilities are from a algorithm such as
#' \code{\link{nbProbabilities}}, then it is recommeneded to use a clustering method
#' and only include the top cluster of infectors for infectees which have such a cluster.
#' This can be specified by using the \code{clustMethod} and \code{cutoff} arguments which
#' are passed into \code{\link{clusterInfectors}}. See the details of this function for
#' a description of the different clustering methods.
#' 
#' The method can be performed with any generation/serial interval distribution,
#' but this version of this function assumes that the generation/serial interval has a gamma distribution.
#' The function does allow for a shifted gamma distribution. The \code{shift} argument
#' defines how much the gamma distribution should be shifted. Any observed generation/serial intervals
#' that are less than this shift will have probability 0. This parameter should be used if 
#' there is a clinical lower bound for the possible generation/serial interval. If this argument
#' is not specified then an unmodified gamma function is used. The units of the
#' estimated gamma distribution will be defined by the units of the provided
#' \code{<timeDiffVar>} column. The value of the \code{shift} should be in the same units.
#' 
#' The algorithm requires initial parameters which should be specified as a vector: 
#' \code{c(<shape>, <scale>)}. These parameters should result in a gamma distribution
#' that is on the desired scale, set by the \code{<timeDiffVar>} column.
#' 
#' If \code{bootSamples > 0}, bootstrap confidence intervals will be estimated for
#' both the shape and scale parameters as well as the mean, median, and mode of the
#' distribution using cluster bootstrapping.
#' 
#' 
#' 
#' @param df The name of the dateset with transmission probabilities (column \code{pVar}),
#' individual IDs (columns \code{<indIDVar>.1} and \code{<indIDVar>.2}), and difference
#' in time between the pair of cases (column \code{timeDiffVar})
#' @param indIDVar The name (in quotes) of the individual ID columns
#' (data frame \code{df} must have variables called \code{<indIDVar>.1}
#'  and \code{<indIDVar>.2}).
#' @param timeDiffVar The name (in quotes) of the column with the difference
#' in time between infection (generation interval) or symptom onset (serial interval) for the
#' pair of cases. The units of this variable (hours, days, years) defines the units of the 
#' resulting distribution.
#' @param pVar The column name (in quotes) of the transmission probabilities.
#' @param clustMethod The method used to cluster the infectors; one of 
#' \code{"none", "n", "kd", "hc_absolute", "hc_relative"} where \code{"none"} or
#' not specifying a value means use all pairs with no clustering
#' (see \code{\link{clusterInfectors}} for detials on clustering methods).
#' @param cutoffs A vector of cutoffs for clustering (see \code{\link{clusterInfectors}}).
#' If more than one cutoff is provided, a pooled estimate will also be provided.
#' @param initialPars A vector of length two with the shape and scale 
#' to initialize the gamma distribution parameters.
#' @param shift A value in the same units as \code{timeDiffVar} that the
#' gamma distribution should be shifted. The default value of 0 is an 
#' unmodifed gamma distribution.
#' @param epsilon The difference between successive estimates of the shape and
#' scale parameters that indicates convergence.
#' @param bootSamples The number of bootstrap samples; if 0, then no confidence intervals
#' are calculated.
#' @param alpha The alpha level for the confidence intervals.
#' @param progressBar A logical indicating if a progress bar should be printed (default is TRUE).
#' 
#'
#' @return A data frame with one row and the following columns:
#' \itemize{
#'    \item \code{nIndividuals} - the number of infectees who have intervals included in the estimate.
#'    \item \code{pCluster} - the proportion of cases who have intervals included in the estimate.
#'    \item \code{nInfectors} - the average number of infectors in the top cluster.
#'    \item \code{shape} - the shape of the estimated gamma distribution for the interval
#'    \item \code{scale} - the scale of the estimated gamma distribution for the interval
#'    \item \code{meanSI} - the mean of the estimated gamma distribution for the interval 
#'    (\code{shape * scale + shift})
#'    \item \code{medianSI} - the median of the estimated gamma distribution for the interval
#'    (\code{qgamma(0.5, shape, scale) + shift)})
#'    \item \code{sdSI} - the standard deviation of the estimated gamma distribution for
#'    the interval (\code{shape * scale ^ 2})
#'  }
#'  If bootSamples > 0, then the data frame also includes the following columns:
#'  \itemize{
#'     \item \code{shapeCILB} and \code{shapeCIUB} - lower bound and upper bounds
#'     of the bootstrap confidence interval for the shape parameter
#'     \item \code{scaleCILB} and \code{scaleCIUB} - lower bound and upper bounds
#'     of the bootstrap confidence interval for the scale parameter
#'     \item \code{meanCILB} and \code{meanCIUB} - lower bound and upper bounds
#'     of the bootstrap confidence interval for the mean of the interval distribution
#'     \item \code{medianCILB} and \code{medianCIUB} - lower bound and upper bounds
#'     of the bootstrap confidence interval for the median of the interval distribution
#'     \item \code{sdCILB} and \code{sdCIUB} - lower bound and upper bounds
#'     of the bootstrap confidence interval for the sd of the interval distribution
#'  }
#' 
#' 
#' @examples
#' 
#' ## First, run the algorithm without including time as a covariate.
#' orderedPair <- pairData[pairData$infectionDiffY > 0, ]
#' 
#' ## Create a variable called snpClose that will define probable links
#' # (<3 SNPs) and nonlinks (>12 SNPs) all pairs with between 2-12 SNPs
#' # will not be used to train.
#' orderedPair$snpClose <- ifelse(orderedPair$snpDist < 3, TRUE,
#'                         ifelse(orderedPair$snpDist > 12, FALSE, NA))
#' table(orderedPair$snpClose)
#' 
#' ## Running the algorithm
#' # NOTE should run with nReps > 1.
#' resGen <- nbProbabilities(orderedPair = orderedPair,
#'                             indIDVar = "individualID",
#'                             pairIDVar = "pairID",
#'                             goldStdVar = "snpClose",
#'                             covariates = c("Z1", "Z2", "Z3", "Z4"),
#'                             label = "SNPs", l = 1,
#'                             n = 10, m = 1, nReps = 1)
#'                             
#' ## Merging the probabilities back with the pair-level data
#' nbResultsNoT <- merge(resGen[[1]], orderedPair, by = "pairID", all = TRUE)
#' 
#' ## Estimating the serial interval
#'
#' # Using hierarchical clustering with a 0.05 absolute difference cutoff
#' estimateSI(nbResultsNoT, indIDVar = "individualID",
#'              timeDiffVar = "infectionDiffY", pVar = "pScaled",
#'              clustMethod = "hc_absolute", cutoff = 0.05, initialPars = c(2, 2))
#'              
#' \donttest{
#' # Using all pairs
#' estimateSI(nbResultsNoT, indIDVar = "individualID",
#'               timeDiffVar = "infectionDiffY", pVar = "pScaled",
#'               clustMethod = "none", initialPars = c(2, 2))
#' 
#'
#' # # Using a shifted gamma distribution:
#' # # not allowing serial intervals of less than 3 months (0.25 years)
#' estimateSI(nbResultsNoT, indIDVar = "individualID",
#'               timeDiffVar = "infectionDiffY", pVar = "pScaled",
#'               clustMethod = "hc_absolute", cutoff = 0.05,
#'               initialPars = c(2, 2), shift = 0.25)
#' 
#' 
#' # # Using multiple cutoffs
#' estimateSI(nbResultsNoT, indIDVar = "individualID",
#'               timeDiffVar = "infectionDiffY", pVar = "pScaled",
#'               clustMethod = "hc_absolute", cutoff = c(0.025, 0.05), initialPars = c(2, 2))
#' }
#' 
#' 
#' ## Adding confidence intervals
#' # NOTE should run with bootSamples > 2.
#' estimateSI(nbResultsNoT, indIDVar = "individualID",
#'              timeDiffVar = "infectionDiffY", pVar = "pScaled",
#'              clustMethod = "hc_absolute", cutoff = 0.05,
#'              initialPars = c(2, 2), shift = 0.25, bootSamples = 2)
#' 
#' @seealso \code{\link{nbProbabilities}} \code{\link{clusterInfectors}}
#'  \code{\link{performPEM}}
#' 
#' @references 
#' Hens N, Calatayud L, Kurkela S, Tamme T, Wallinga J. Robust reconstruction and
#' analysis of outbreak data: influenza A (H1N1) v transmission in a school-based 
#' population. \emph{American Journal of Epidemiology}. 2012 Jul 12;176(3):196-203.
#' 
#' @export
 

estimateSI <- function(df, indIDVar, timeDiffVar, pVar,
                       clustMethod = c("none", "n", "kd", "hc_absolute", "hc_relative"),
                       cutoffs = NULL, initialPars, shift = 0, epsilon = 0.00001,
                       bootSamples = 0, alpha = 0.05, progressBar = TRUE){
  
  #Creating variables with the individual ID
  indIDVar1 <- paste0(indIDVar, ".1")
  indIDVar2 <- paste0(indIDVar, ".2")
  
  #Checking that the named variables are in the data frame
  if(!indIDVar1 %in% names(df)){
    stop(paste0(indIDVar1, " is not in the data frame."))
  }
  if(!indIDVar2 %in% names(df)){
    stop(paste0(indIDVar2, " is not in the data frame."))
  }
  if(!timeDiffVar %in% names(df)){
    stop(paste0(timeDiffVar, " is not in the data frame."))
  }
  if(!pVar %in% names(df)){
    stop(paste0(pVar, " is not in the data frame."))
  }
  
  #If clustMethod is not specified, setting it to "none"
  if(length(clustMethod) > 1){
    clustMethod <- "none"
    warning("No clustMethod was provided so it was set to 'none'")
  }
  if(!clustMethod %in% c("none", "n", "kd", "hc_absolute", "hc_relative")){
    stop("clustMethod must be one of: none, n, kd, hc_absolute, hc_relative")
  }
  
  #Make sure cutoff is provided if clustMethod is not "none"
  if(clustMethod != "none" & is.null(cutoffs)){
    stop("Please provide one or more cutoff values")
  }
  
  df <- as.data.frame(df[, c(indIDVar1, indIDVar2, timeDiffVar, pVar)],
                      stringsAsFactors = FALSE)

  
  #Storing a list of clustering results for all cutoffs
  if(clustMethod != "none"){
    
    clustList <- NULL
    
    if(progressBar == TRUE){
      pb <- utils::txtProgressBar(min = 0, max = length(cutoffs), style = 3) 
    }
    
    for(i in 1:length(cutoffs)){
      cutoff <- cutoffs[i]
      clustList[[i]] <- clusterInfectors(df, indIDVar = indIDVar, pVar = pVar,
                                         clustMethod = clustMethod, cutoff = cutoff)
      
      if(progressBar == TRUE){
        utils::setTxtProgressBar(pb, i) 
      }
    }
  }else{
    df$cluster <- 1
    clustList <- list(df)
    cutoffs <- "none"
  }

  
  #Finding the point estimate for the serial interval parameters
  siEst <- estimateSIPars(clustList, indIDVar = indIDVar, timeDiffVar = timeDiffVar,
                          pVar = pVar, clustMethod = clustMethod, cutoffs = cutoffs,
                          initialPars = initialPars, shift = shift, epsilon = epsilon) 
  
  #If bootSamples > 0, estimate bootstrap confidence intervals on the parameters
  if(bootSamples > 0){
    
    bootSI <- NULL
    
    if(progressBar == TRUE){
      pb2 <- utils::txtProgressBar(min = 0, max = bootSamples, style = 3) 
    }
    
    for(i in 1:bootSamples){
      
      #Find a random sample of infectees
      ids <- base::sample(unique(df[, indIDVar2]),
                    size = length(unique(df[, indIDVar2])),
                    replace = TRUE)
      
      #Adding new ID variable
      idDf <- cbind.data.frame(ids, newID = 1:length(ids))
      names(idDf) <- c(indIDVar2, "newID")
      
      #Function to subset the dataset to just the cases in the bootstrap sample
      bootSample <- function(clustdf){
        
        #Using the fact that joining with duplicates, duplicates the column to
        #replicate and add new ID in one step.
        probsBoot <- dplyr::right_join(clustdf, idDf, by = indIDVar2)
        
        #Renaming the new ID variable to the old ID name so that the function works
        names(probsBoot)[names(probsBoot) == indIDVar2] <- "oldID2" 
        names(probsBoot)[names(probsBoot) == "newID"] <- indIDVar2
        return(probsBoot)
      }

      #Adding all infectors for each case in the bootstrap sample
      probsBootList <- lapply(clustList, bootSample)
      
      siNew <- estimateSIPars(probsBootList, indIDVar = indIDVar, timeDiffVar = timeDiffVar,
                              pVar = pVar, clustMethod = clustMethod, cutoffs = cutoffs,
                              initialPars = initialPars, shift = shift, epsilon = epsilon) 
      
      bootSI <- dplyr::bind_rows(bootSI, siNew)
      
      if(progressBar == TRUE){
        utils::setTxtProgressBar(pb2, i)
      }

    }
    
    #Finding the CI bounds
    siQuantsL <- by(bootSI,
                   INDICES = list(bootSI$cutoff),
                   FUN = function(x){
                     data.frame("cutoff" = unique(x$cutoff),
                                "shapelb" = stats::quantile(x$shape, 1-alpha/2, na.rm = TRUE),
                                "shapeub" = stats::quantile(x$shape, alpha/2, na.rm = TRUE),
                                "scalelb" = stats::quantile(x$scale, 1-alpha/2, na.rm = TRUE),
                                "scaleub" = stats::quantile(x$scale, alpha/2, na.rm = TRUE),
                                "meanlb" = stats::quantile(x$meanSI, 1-alpha/2, na.rm = TRUE),
                                "meanub" = stats::quantile(x$meanSI, alpha/2, na.rm = TRUE),
                                "medianlb" = stats::quantile(x$medianSI, 1-alpha/2, na.rm = TRUE),
                                "medianub" = stats::quantile(x$medianSI, alpha/2, na.rm = TRUE),
                                "sdlb" = stats::quantile(x$sdSI, 1-alpha/2, na.rm = TRUE),
                                "sdub" = stats::quantile(x$sdSI, alpha/2, na.rm = TRUE))
                   })
    #Suppressing character to factor warnings
    suppressWarnings(siQuants <- do.call(dplyr::bind_rows, siQuantsL))
    siQuants$cutoff <- as.character(siQuants$cutoff)
    
    siCI <- dplyr::full_join(siEst, siQuants, by = "cutoff")
    siCI$shapeCILB <- siCI$shape - (siCI$shapelb - siCI$shape)
    siCI$shapeCIUB <- siCI$shape - (siCI$shapeub - siCI$shape)
    siCI$scaleCILB <- siCI$scale - (siCI$scalelb - siCI$scale)
    siCI$scaleCIUB <- siCI$scale - (siCI$scaleub - siCI$scale)
    siCI$meanCILB <- siCI$meanSI - (siCI$meanlb - siCI$meanSI)
    siCI$meanCIUB <- siCI$meanSI - (siCI$meanub - siCI$meanSI)
    siCI$medianCILB <- siCI$medianSI - (siCI$medianlb - siCI$medianSI)
    siCI$medianCIUB <- siCI$medianSI - (siCI$medianub - siCI$medianSI)
    siCI$sdCILB <- siCI$sdSI - (siCI$sdlb - siCI$sdSI)
    siCI$sdCIUB <- siCI$sdSI - (siCI$sdub - siCI$sdSI)
    
    siCI <- siCI[, !names(siCI) %in% c("shapelb", "shapeub", "scalelb", "scaleub", "meanlb",
                                       "meanub", "medianlb", "medianub", "sdlb", "sdub")]
    
    return(siCI)
    
  }else{
    return(siEst)
  }
}
  




# Function run inside estimateSI to allow for bootstrap CI; not run on its own.
estimateSIPars <- function(clustL, indIDVar, timeDiffVar, pVar, clustMethod,
                       cutoffs, initialPars, shift, epsilon){
  
  #Creating variables with the individual ID
  indIDVar1 <- paste0(indIDVar, ".1")
  indIDVar2 <- paste0(indIDVar, ".2")
  
  #Estimating SI for all given cutoffs
  siData <- NULL
  for(i in 1:length(clustL)){
    
    cutoff <- cutoffs[i]
    clustRes <- clustL[[i]]
    
    #Restricting to just the top cluster to be used for estimation
    topClust <- clustRes[clustRes$cluster == 1, ]
    
    #Finding the number of infectees with a top cluster used to estimate the serial interval
    nIndividuals <- length(unique(topClust[, indIDVar2]))
    #Finding the proportion of cases with a top cluster
    pCluster <- nIndividuals / length(unique(clustRes[, indIDVar2]))
    #Finding the average number of infectors per infectee in the top cluster
    nInfectors <- mean(as.numeric(table(topClust[, indIDVar2])))
    
    
    #Estimating the serial interval parameters using the PEM algorithm
    pars <- NULL
    if(nIndividuals >= 10){
      pars <- performPEM(df = topClust, indIDVar = indIDVar,
                         timeDiffVar = timeDiffVar,
                         pVar = pVar, initialPars = initialPars,
                         shift = shift, epsilon = epsilon)
      
    }else{
      message(paste0("With ", clustMethod, " and ", cutoff,
                   ", fewer than 10 individuals would be used for estimation so this cuttof is excluded"))
      pars <- cbind.data.frame("shape" = NA, "scale" = NA)
    }
    
    siDataTemp <- cbind.data.frame(clustMethod, cutoff = as.character(cutoff),
                                   nIndividuals, nInfectors, pCluster, pars,
                                   stringsAsFactors = FALSE)
    siData <- dplyr::bind_rows(siData, siDataTemp)
  }
  
  #Adding summary measures
  siData$meanSI <- siData$shape * siData$scale + shift
  siData$medianSI <- stats::qgamma(0.5, shape = siData$shape,
                                   scale = siData$scale) + shift
  siData$sdSI <- sqrt(siData$shape * siData$scale ^ 2)
  
  #Adding pooled estimates if more than one cutoff is provided
  if(length(cutoffs) > 1){
    siDataPooled <- as.data.frame(t(sapply(siData[, c("meanSI", "medianSI", "sdSI")],
                                           mean, na.rm = TRUE)), stringsAsFactors = FALSE)
    siDataPooled$clustMethod <- clustMethod
    siDataPooled$cutoff <- "pooled"
    siData <- dplyr::bind_rows(siData, siDataPooled) 
  }
    
  return(siData)
}




#' Executes the PEM algorthim to estimate the generation/serial interval distribution
#'
#' The function \code{performPEM} uses relative transmission probabilities to estimate
#' the generation/serial interval distribution
#' 
#' This function is meant to be called by \code{\link{estimateSI}}
#' which estimates the generation/serial interval distribution as well as clustering the
#' probabilities, but can be called directly. The main reason to call \code{performPEM}
#' directly is for the \code{plot} argument. Setting this argument to \code{TRUE}
#' will produce a plot of the shape and scale parameters at each iteration.
#' For more details on the PEM algorithm see \code{\link{estimateSI}}.
#' 
#' 
#' @param df The name of the dateset with transmission probabilities.
#' @param indIDVar The name (in quotes) of the individual ID columns
#' (data frame \code{df} must have variables called \code{<indIDVar>.1}
#'  and \code{<indIDVar>.2}).
#' @param timeDiffVar The name (in quotes) of the column with the difference
#' in time between infection (generation interval) or symptom onset (serial interval) for the
#' pair of cases. The units of this variable (hours, days, years) defines the units of the 
#' resulting distribution.
#' @param pVar The column name (in quotes) of the transmission probabilities.
#' @param initialPars A vector of length two with the shape and scale 
#' to initialize the gamma distribution parameters.
#' @param shift A value in the same units as \code{timeDiffVar} that the
#' gamma distribution should be shifted. The Default value of 0 is an 
#' unmodifed gamma distribution.
#' @param epsilon The difference between successive estimates of the shape and
#' scale parameters that indicates convergence.
#' @param plot A logical indicating if a plot should be printed showing the
#' parameter estimates at each iteration.
#' 
#'
#' @return A data frame with one row and the following columns:
#' \itemize{
#'    \item \code{nIndividuals} - the number of infectees who have intervals included
#'    in the SI estimate.
#'    \item \code{shape} - the shape of the estimated gamma distribution for the interval.
#'    \item \code{scale} - the scale of the estimated gamma distribution for the interval.
#'    \item \code{meanSI} - the mean of the estimated gamma distribution for the interval 
#'    (\code{shape * scale + shift}).
#'    \item \code{medianSI} - the median of the estimated gamma distribution for the interval
#'    (\code{qgamma(0.5, shape, scale) + shift)}).
#'    \item \code{sdSI} - the standard deviation of the estimated gamma distribution for
#'    the interval (\code{shape * scale ^ 2})
#'  }
#' 
#' 
#' @examples
#' 
#' ## First, run the algorithm without including time as a covariate.
#' orderedPair <- pairData[pairData$infectionDiffY > 0, ]
#' 
#' ## Create a variable called snpClose that will define probable links
#' # (<3 SNPs) and nonlinks (>12 SNPs) all pairs with between 2-12 SNPs
#' # will not be used to train.
#' orderedPair$snpClose <- ifelse(orderedPair$snpDist < 3, TRUE,
#'                         ifelse(orderedPair$snpDist > 12, FALSE, NA))
#' table(orderedPair$snpClose)
#' 
#' ## Running the algorithm
#' #NOTE should run with nReps > 1.
#' resGen <- nbProbabilities(orderedPair = orderedPair,
#'                             indIDVar = "individualID",
#'                             pairIDVar = "pairID",
#'                             goldStdVar = "snpClose",
#'                             covariates = c("Z1", "Z2", "Z3", "Z4"),
#'                             label = "SNPs", l = 1,
#'                             n = 10, m = 1, nReps = 1)
#'                             
#' ## Merging the probabilities back with the pair-level data
#' nbResultsNoT <- merge(resGen[[1]], orderedPair, by = "pairID", all = TRUE)
#' 
#' ## Estimating the serial interval
#' 
#' \donttest{
#' # Using all pairs and plotting the parameters
#'   performPEM(nbResultsNoT, indIDVar = "individualID", timeDiffVar = "infectionDiffY",
#'   pVar = "pScaled", initialPars = c(2, 2), shift = 0, plot = TRUE)
#' }
#'
#'
#' # Clustering the probabilities first
#' allClust <- clusterInfectors(nbResultsNoT, indIDVar = "individualID", pVar = "pScaled",
#'                             clustMethod = "hc_absolute", cutoff = 0.05)
#' 
#' performPEM(allClust[allClust$cluster == 1, ], indIDVar = "individualID",
#'            timeDiffVar = "infectionDiffY", pVar = "pScaled",
#'            initialPars = c(2, 2), shift = 0, plot = TRUE)
#'
#' \donttest{
#' # The above is equivalent to the following code using the function estimateSI()
#' # though the plot will not be printed and more details will be added
#' estimateSI(nbResultsNoT, indIDVar = "individualID", timeDiffVar = "infectionDiffY",
#'            pVar = "pScaled", clustMethod = "hc_absolute", cutoff = 0.05,
#'            initialPars = c(2, 2))
#' }
#' 
#' 
#' @seealso \code{\link{nbProbabilities}} \code{\link{clusterInfectors}}
#'  \code{\link{performPEM}}
#' 
#' @references 
#' Hens N, Calatayud L, Kurkela S, Tamme T, Wallinga J. Robust reconstruction and
#' analysis of outbreak data: influenza A (H1N1) v transmission in a school-based 
#' population. \emph{American Journal of Epidemiology}. 2012 Jul 12;176(3):196-203.
#' 
#' @export


performPEM <- function(df, indIDVar, timeDiffVar, pVar, initialPars,
                       shift = 0, epsilon = 0.00001, plot = FALSE){
  
  df <- as.data.frame(df)
  #Creating variables with the individual ID
  indIDVar1 <- paste0(indIDVar, ".1")
  indIDVar2 <- paste0(indIDVar, ".2")
  
  #Naming the columns correctly
  dfNew <- as.data.frame(df[, c(indIDVar1, indIDVar2, timeDiffVar, pVar)])
  
  #Initializing the data frame with the starting parameters
  i <- 1
  params <- as.data.frame(t(c(i, initialPars[1], initialPars[2], NA, NA)))
  names(params) <- c("iteration", "shape", "scale", "diff.shape", "diff.scale")
  
  #Continue to update the parameters until the difference between the current parameters
  #and the one's in the previous iteration is less than epsilon
  while(is.na(params[i, "diff.shape"]) | is.na(params[i, "diff.scale"]) |
        params[i, "diff.shape"] > epsilon | params[i, "diff.scale"] > epsilon){
    
    #Updating the parameter estimate
    paramsNew <- updateParams(df = dfNew, indIDVar = indIDVar, timeDiffVar = timeDiffVar,
                              pVar = pVar, currentPars = c(params[i, "shape"],
                                                           params[i, "scale"]),
                              shift = shift)
    
    #Adding the new estimate to the data frame
    params <- dplyr::bind_rows(params, paramsNew)
    params[i+1, "iteration"] <- i+1
    params[i+1, "diff.shape"] <- abs(params[i+1, "shape"] - params[i, "shape"])
    params[i+1, "diff.scale"] <- abs(params[i+1, "scale"] - params[i, "scale"])
    i <- i + 1

    if(i == 100){
      warning(paste0("Results are for 100 iterations, which did not meet epsilon difference \n",
                     "Shape difference was ", signif(params[i, "diff.shape"], 2),
                     " and scale difference was ", signif(params[i, "diff.scale"], 2)))
      break()
      }
  }
  paramsFinal <- params[nrow(params), c("shape", "scale")]
  
  #If asked for, plot the parameters over the optimization
  if(plot == TRUE){
    paramsLong <- params[, c("shape", "scale", "iteration")]
    paramsLong <- tidyr::gather(paramsLong, !!rlang::sym("parameter"),
                                !!rlang::sym("value"), -!!rlang::sym("iteration"))
    p <- ggplot2::ggplot(data = paramsLong) +
      ggplot2::geom_point(ggplot2::aes(x = !!rlang::sym("iteration"), 
                                       y = !!rlang::sym("value"),
                                       color = !!rlang::sym("parameter")))
    print(p)
  }
  
  return(paramsFinal)
}


#Function to update the gamma parameters using the current parameters
updateParams <- function(df, indIDVar, timeDiffVar, pVar, currentPars, shift = 0){
  
  df <- as.data.frame(df)
  #Creating variables with the individual ID
  indIDVar1 <- paste0(indIDVar, ".1")
  indIDVar2 <- paste0(indIDVar, ".2")
  
  #Calculating serial interval distribution times the probability
  df$gp <- dgammaS(x = df[, timeDiffVar], shift = shift,
                   shape = currentPars[1], scale = currentPars[2]) * df[, pVar]
  
  #Calculating the total of this value for each person
  total.gp <- stats::aggregate(df$gp, by = list(df[, indIDVar2]),
                           sum, na.rm = TRUE)
  names(total.gp) <- c(indIDVar2, "gpTotal")

  #Adding back this total to get the new scaled probability (pij)
  df2 <- merge(df, total.gp, by = indIDVar2, all = TRUE)
  df2$pij <- ifelse(df2$gpTotal != 0, df2$gp / df2$gpTotal, 0)
  
  si <- stats::optim(par = c(currentPars[1], currentPars[2]),
                     logl_pem, df = df2, timeDiffVar = timeDiffVar,
                     shift = shift, method = "L-BFGS-B",
                     lower = c(.0005, 0.0005), upper = c(Inf, Inf))
  names(si$par) <- c("shape", "scale")
  
  return(si$par)
}


#Function to calculate the density for the shifted gamma
dgammaS <- function(x, shift, shape, scale){
  shape <- as.numeric(shape)
  scale <- as.numeric(scale)
  d <- ifelse(x - shift <= 0, 0, stats::dgamma(x - shift, shape = shape, scale = scale))
  return(d)
}

#Log likelihood function for PEM algorithm for SI parameters
logl_pem <- function(pars, df, timeDiffVar, shift = 0){
  df$g <- dgammaS(x = df[, timeDiffVar], shape = pars[1], scale = pars[2],
                  shift = shift)
  df$logg <- ifelse(df$g == 0, 0, log(df$g))
  -sum(df$pij * df$logg, na.rm = TRUE)
}

#Log likelihood function for basic MLE for SI parameters
logl_si <- function(params, df, timeDiffVar, shift = 0){
  df$g <- dgammaS(x = df[, timeDiffVar], shape = params[1], scale = params[2],
                  shift = shift)
  df$logg <- ifelse(df$g == 0, 0, log(df$g))
  -sum(df$logg, na.rm = TRUE)
}
