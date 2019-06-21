
#' Perform Naive Bayes
#'
#' Calculates predicted probabilities using naive Bayes
#'
#' @param training training dataset
#' @param validation validation dataset
#' @param covariates vector of the covariate variable names
#' @param goldStdVar the variable name (in quotes) that will define linking status
#' @param weighting should deep frequency covariate weighting be applied
#' @param label optional label for the run (default is NULL)
#'
#' @return List containing two dataframes: "probs" with pairdata with an extra column with the
#'    probabilities and "coeff" with the coefficient values.
#'
#' @examples
#' #Insert example here
#'
#' @export


performNB <- function(training, validation, covariates,
                      goldStdVar, weighting=FALSE, label){
  
  #Create label
  if (weighting == FALSE){labelFull <- paste0(label, goldStdVar, "NB")}
  if (weighting == TRUE){labelFull <- paste0(label, goldStdVar, "NBDFW")}
  

  #### Function for CFS Search Algorithms ####
  
  calcEntropy <- function(subset){
    
    k <-  length(subset)
    
    if(k > 1){
      #Finding all pairs of variables
      pairs <- as.data.frame(t(combn(subset, 2)), stringsAsFactors = FALSE)
      names(pairs) <- c("Var1", "Var2")
      pairs <- pairs %>% mutate(comb = paste(Var1, Var2, sep = "~"))
      
      #Finding symmertric uncertainty (adjusted information gain) for all pairs of features
      Rff <- map_df(pairs$comb, symmetrical.uncertainty, data=training) 
      #Mean of all feature-feature correlations
      rff <- mean(Rff[,1])
    }else{rff <- 0}
    
    model <- as.simple.formula(subset, "linked")
    Rcf <- symmetrical.uncertainty(model, training)
    #Mean of all classification-feature correlations
    rcf <- mean(Rcf[,1])
    
    #Calculating merit of the subset
    merit <- k * rcf / sqrt(k + k * (k - 1) * rff)
    #print(subset)
    #print(merit)
    return(merit)
  }

    
  #### Naive Bayes Method ####
  
    #Setting up weight table
    varTable <- as.data.frame(covariates)
    varTable <- varTable %>% mutate(variable = as.character(covariates))
    
    if(weighting == TRUE){
      #Running CFS with with two different methods for calculating correlation
      model <- best.first.search(covariates, calcEntropy)
      #Setting weights to 2 if the variable is chosen and 1 if it is not chosen
      varTable <- varTable %>% mutate(weight = ifelse(variable %in% model, 2, 1))
    }else{varTable$weight <- 1}
    
    
    #Creating the results dataframe
    results <- validation
    #Finding proportion of links/non-links
    classTab <- prop.table(table(training[, "linked"]) + 1)
    coeff <- NULL
    
    #Looping through all covariates
    for(i in 1:length(covariates)){
      
      #Extracting covariate name
      Var <- covariates[i]
      #Determining the weight for that covariate
      W <- varTable %>% filter(variable == Var) %>% pull(weight)
      #Creating a table with proportions in each level of the covariate from training data
      Tab <- prop.table(W * table(training[, Var], training[, "linked"]) + 1, 2)
      
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
      coeff$label <- labelFull
    }
    
    #Calculating numerator and denominator for the probability calculation
    results$link <- apply(results[, grepl("_T", names(results))], 1, prod) * classTab[2] 
    results$nonlink <- apply(results[, grepl("_F", names(results))], 1, prod) * classTab[1]
    
    #Calculating probability of link
    probs <- (results
                %>% mutate(p = link / (link + nonlink))
                %>% bind_rows(training)
                %>% mutate(nLinksTrain = sum(training$linked == TRUE),
                           label = labelFull)
                %>% select(label, edgeID, p, nLinksTrain)
    )
    
    return(list(probs, coeff))
}