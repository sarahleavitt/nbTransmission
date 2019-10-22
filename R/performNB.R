
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
                      covariates, l = 1){
  
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
    varTable$weight <- 1
    
    #Creating the results dataframe
    results <- validation
    #Finding proportion of links/non-links
    classTab <- prop.table(table(training[, goldStdVar]) + l)
    coeff <- NULL
    
    #Looping through all covariates
    for(i in 1:length(covariates)){
      
      #Extracting covariate name
      Var <- covariates[i]
      #Determining the weight for that covariate
      W <- varTable %>% filter(variable == Var) %>% pull(weight)
      #Creating a table with proportions in each level of the covariate from training data
      Tab <- prop.table(W * table(training[, Var], training[, goldStdVar]) + l, 2)
      
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
