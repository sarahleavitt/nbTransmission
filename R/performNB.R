
#' Perform Naive Bayes
#'
#' Calculates predicted probabilities using naive Bayes
#' 
#' Add details section
#'
#' @param training The training dataset name.
#' @param validation The validation dataset name.
#' @param obsIDVar The variable name (in quotes) of the observation ID variable.
#' @param goldStdVar The variable name (in quotes) that will define event status
#'  in the training dataset.
#' @param covariates A character vector containing the covariate variable names.
#' @param l Laplace smoothing parameter that is added to each cell (default is 1).
#'
#' @return List containing two dataframes: "probs" containing both training and validation 
#' observations with the two columns: obsIDVar and p and "coeff" with the coefficient values.
#'
#' @examples
#' #Insert example here
#' 
#' @export



performNB <- function(training, validation, obsIDVar, goldStdVar, 
                      covariates, l = 1){
  
  #### Checking variable names ####
  
  #Checking that the named variables are in the dataframe
  if(!obsIDVar %in% names(training) | !obsIDVar %in% names(validation)){
    stop(paste0(obsIDVar, " is not in the dataframe."))
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
  
  
  #Making sure there are both events and non-event observations.
  #If not return NA for probabilities and print a warning
  if(sum(training[, goldStdVar] == TRUE, na.rm = TRUE) == 0 |
     sum(training[, goldStdVar] == FALSE, na.rm = TRUE) == 0){
    
    #Setting probability of a event to 0
    probs <- validation
    probs$p <- NA
    if(!"p" %in% names(training)){training$p <- NA}
    probs <- rbind(probs[, c("p", obsIDVar)], training[, c("p", obsIDVar)])
    coeff <- NULL
    
    warning("No events or no non-events in training set")
    
  }else{
    
    #Setting up weight table
    #(currently irrelevant but would allow for weighting in subsequent versions)
    varTable <- as.data.frame(covariates)
    varTable$variable <- as.character(covariates)
    varTable$weight <- 1
    
    #Creating the results dataframe which is a copy of the validation dataframe
    results <- validation
    #Finding proportion of events/non-events in the trianing dataset
    classTab <- prop.table(table(training[, goldStdVar]) + l)
    
    #Initializing the coefficient dataframe
    coeff <- data.frame("level" = character(), "ratio" = numeric())
    
    #Looping through all covariates and finding frequencies in training data
    for(i in 1:length(covariates)){
      
      #Extracting covariate name
      Var <- covariates[i]
      #Determining the weight for that covariate
      W <- varTable[varTable$variable == Var, "weight"]
      #Creating a table with proportions in each level of the covariate from training data
      Tab <- prop.table(W * table(training[, Var], training[, goldStdVar]) + l, 2)
      
      #Calculating P(X|Outcome = TRUE, W) and P(X|Outcome = FALSE, W)
      #First creating a variable with the value of this covariate
      results$covariate <- results[, Var]
      #1 if missing, otherwise it is the proportion in the second column and the row
      #associated with the factor level of the covariate value.
      results[, paste0(Var, "_T")] <- ifelse(is.na(results$covariate), 1,
                                             Tab[as.numeric(results$covariate), 2]) ^ W
      results[, paste0(Var, "_F")] <- ifelse(is.na(results$covariate), 1,
                                             Tab[as.numeric(results$covariate), 1]) ^ W
      
      #Saving the P(X|Outcome = TRUE, W) / P(X|Outcome = FALSE, W)
      ratio <- Tab[, 2] / Tab[, 1]
      level <- paste(Var, names(ratio), sep = ":")
      cTemp <- cbind.data.frame(level, ratio, stringsAsFactors = FALSE)
      coeff <- rbind(coeff, cTemp)
    }
    
    #Calculating numerator and denominator for the probability calculation
    results$event <- apply(results[, grepl("_T", names(results))], 1, prod) * classTab[2] 
    results$nonevent <- apply(results[, grepl("_F", names(results))], 1, prod) * classTab[1]
    
    #Calculating probability of event
    probs <- results
    probs$p <- probs$event / (probs$event + probs$nonevent)
    
    #Combining training and validation datasets
    if(!"p" %in% names(training)){training$p <- NA}
    probs <- rbind(probs[, c("p", obsIDVar)], training[, c("p", obsIDVar)])
  }
  
  return(list(probs, coeff))
}
