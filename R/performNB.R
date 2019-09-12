
#' Perform Naive Bayes
#'
#' Calculates predicted probabilities using naive Bayes
#' 
#' Add details section
#'
#' @param training The training dataset name.
#' @param validation The validation dataset name.
#' @param edgeIDVar The variable name (in quotes) of the edge ID variable.
#' @param goldStdVar The variable name (in quotes) that will define linking status.
#' @param covariates A character vector containing the covariate variable names.
#' @param l Laplace smoothing parameter that is added to each cell (default is 1).
#' @param nbWeighting A logical scalar. Do you want to use deep frequency weighting in NB (default is FALSE)?
#'
#' @return List containing two dataframes: "probs" with pairdata with an extra column with the
#'    probabilities and "coeff" with the coefficient values.
#'
#' @examples
#' #Insert example here
#'
#' @import dplyr
#' 
#' @export


performNB <- function(training, validation, edgeIDVar, goldStdVar, 
                      covariates, l = 1, nbWeighting=FALSE){
  
  #### Checking variable names ####
  
  #Checking that the named variables are in the dataframe
  if(!edgeIDVar %in% names(training) | !edgeIDVar %in% names(validation)){
    stop(paste0(edgeIDVar, " is not in the dataframe."))
  }
  if(!goldStdVar %in% names(training)){
    stop(paste0(goldStdVar, " is not in the dataframe."))
  }
  #Checking that the covariates are in the dataframe
  covarTestT <- c(covariates %in% names(training), covariates %in% names(validation))
  if(FALSE %in% covarTestT){
    stop("At least one of the covariates is not in the input dataframes.")
  }

  #Checking that all of the covariates are factors
  covarDf <- training[, covariates]
  notFactor <- names(covarDf)[!sapply(covarDf, is.factor)]
  notFactorC <- paste0(notFactor, collapse = ", ")
  if(FALSE %in% sapply(covarDf, is.factor)){
    stop(paste0(notFactorC, " are not factors in the training data"))
  }
  covarDf <- validation[, covariates]
  notFactor <- names(covarDf)[!sapply(covarDf, is.factor)]
  notFactorC <- paste0(notFactor, collapse = ", ")
  if(FALSE %in% sapply(covarDf, is.factor)){
    stop(paste0(notFactorC, " are not factors in the training data"))
  }
  
  
  #Making sure there are both linked and nonlinked pairs.
  #If not return NA for probabilities and print a warning
  if(sum(training[, goldStdVar] == TRUE, na.rm = TRUE) == 0 |
     sum(training[, goldStdVar] == FALSE, na.rm = TRUE) == 0){
    
    #Calculating probability of link
    probs <- (validation
              %>% mutate(p = NA)
              %>% bind_rows(training)
    )
    probs <- probs[, c("p", edgeIDVar)]
    coeff <- NULL
    
    warning("No events or non-events in training set")
    
  }else{
    
    #Setting up weight table
    varTable <- as.data.frame(covariates)
    varTable <- varTable %>% mutate(variable = as.character(covariates))
    
    if(nbWeighting == TRUE){
      #Running CFS with with two different methods for calculating correlation
      model <- FSelector::best.first.search(covariates, calcEntropy)
      #Setting weights to 2 if the variable is chosen and 1 if it is not chosen
      varTable <- varTable %>% mutate(weight = ifelse(variable %in% model, 2, 1))
    }else{varTable$weight <- 1}
    
    
    #Creating the results dataframe
    results <- validation
    #Finding proportion of links/non-links
    classTab <- prop.table(table(training[, "linked"]) + l)
    coeff <- NULL
    
    #Looping through all covariates
    for(i in 1:length(covariates)){
      
      #Extracting covariate name
      Var <- covariates[i]
      #Determining the weight for that covariate
      W <- varTable %>% filter(variable == Var) %>% pull(weight)
      #Creating a table with proportions in each level of the covariate from training data
      Tab <- prop.table(W * table(training[, Var], training[, "linked"]) + l, 2)
      
      #Determining P(Y|L = TRUE, W) and P(Y|L = FALSE, W)
      #First creating a variable with the value of this covariate
      results$covariate <- results[, Var]
      #1 if missing, otherwise it is the proportion in the second column and the row
      #associated with the factor level of the covariate value.
      results[, paste0(Var, "_T")] <- ifelse(is.na(results$covariate), 1,
                                             Tab[as.numeric(results$covariate), 2]) ^ W
      results[, paste0(Var, "_F")] <- ifelse(is.na(results$covariate), 1,
                                             Tab[as.numeric(results$covariate), 1]) ^ W
      
      #Saving the P(Y|L = TRUE) / P(Y|L = FALSE)
      ratio <- Tab[, 2] / Tab[, 1]
      level <- paste(Var, names(ratio), sep = ":")
      cTemp <- cbind.data.frame(level, ratio, stringsAsFactors = FALSE)
      coeff <- bind_rows(coeff, cTemp)
    }
    
    #Calculating numerator and denominator for the probability calculation
    results$link <- apply(results[, grepl("_T", names(results))], 1, prod) * classTab[2] 
    results$nonlink <- apply(results[, grepl("_F", names(results))], 1, prod) * classTab[1]
    
    #Calculating probability of link
    probs <- (results
              %>% mutate(p = link / (link + nonlink))
              %>% bind_rows(training)
    )
    probs <- probs[, c("p", edgeIDVar)]
  }
  
  return(list(probs, coeff))
}


#### Function for CFS Search Algorithms ####

calcEntropy <- function(subset){
  
  k <-  length(subset)
  
  if(k > 1){
    #Finding all pairs of variables
    pairs <- as.data.frame(t(utils::combn(subset, 2)), stringsAsFactors = FALSE)
    names(pairs) <- c("Var1", "Var2")
    pairs <- pairs %>% tidyr::unite(comb, Var1, Var2, sep = "~", remove = FALSE)
    
    #Finding symmertric uncertainty (adjusted information gain) for all pairs of features
    Rff <- purrr::map_df(pairs$comb, FSelector::symmetrical.uncertainty, data=training) 
    #Mean of all feature-feature correlations
    rff <- mean(Rff[,1])
  }else{rff <- 0}
  
  model <- FSelector::as.simple.formula(subset, "linked")
  Rcf <- FSelector::symmetrical.uncertainty(model, training)
  #Mean of all classification-feature correlations
  rcf <- mean(Rcf[,1])
  
  #Calculating merit of the subset
  merit <- k * rcf / sqrt(k + k * (k - 1) * rff)
  
  return(merit)
}

