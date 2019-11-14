
#' Perform Naive Bayes
#'
#' The function \code{performNB} Calculates the posterior probabilities of a dichotomous class
#' variable given a set of covariates using Bayes rule.
#' 
#' The main purpose of this function is to be used by \code{\link{nbProbabilities}} to 
#' estimate the relative transmission probability between individuals in an infectious
#' disease outbreak. However, it can be used more generally to estimate the probability
#' of any dichotomous outcome given a set of categorical covariates.
#' 
#' The function needs a training dataset with the outcome variable (\code{goldStdVar})
#' which is \code{TRUE} for those who have the value of interest and \code{FALSE}
#' for those who do not. The probability of having the outcome 
#' (\code{<goldStdVar> = TRUE}) is predicted in the prediction dataset.
#' 
#'
#' @param training The training dataset name.
#' @param prediction The prediction dataset name.
#' @param obsIDVar The variable name (in quotes) of the observation ID variable.
#' @param goldStdVar The variable name (in quotes) of the outcome in the training dataset
#' (needs to be a logical variable with value \code{TRUE} for observations with
#'  the outcome of interest.
#' @param covariates A character vector containing the covariate variable names.
#' All covariates need to be categorical factor variables.
#' @param l Laplace smoothing parameter that is added to each cell
#' (a value of 0 indicates no smoothing).
#'
#' @return List containing two dataframes: 
#' \enumerate{
#'   \item \code{probabilities} - a dataframe combining \code{training} and \code{prediction}
#'    with predictied probabilities for the \code{prediction} dataframe. Column names:
#'      \itemize{
#'        \item \code{<obsIDVar>} - the observation ID with the name specified
#'        \item \code{p} - the probability that \code{<goldStdVar> = TRUE} for observations in the
#'        \code{prediction} dataset.
#'      }
#'   \item \code{estimates} - a dataframe with the effect estimates derived from the training dataset.
#'   Column names:
#'      \itemize{
#'        \item \code{level} - the covariate name and level
#'        \item \code{odds} - the value of the likelihood odds (will change to OR)
#'      }
#' }
#' @seealso \code{\link{nbProbabilities}}
#'
#' @examples
#' ## Use iris dataset and predict if a flower is of the specices "virginica".
#' 
#' data(iris)
#' irisNew <- iris
#' ## Creating an id variable
#' irisNew$id <- seq(1:nrow(irisNew))
#' ## Creating logical variable indicating if the flower is of the species virginica
#' irisNew$spVirginica <- irisNew$Species == "virginica"
#'
#' ## Creating categorical/factor versions of the covariates
#' irisNew$Sepal.Length.Cat <- factor(cut(irisNew$Sepal.Length, c(0, 5, 6, 7, Inf)),
#'                                  labels = c("<=5.0", "5.1-6.0", "6.1-7.0", "7.1+"))
#'
#' irisNew$Sepal.Width.Cat <- factor(cut(irisNew$Sepal.Width, c(0, 2.5, 3, 3.5, Inf)),
#'                                  labels = c("<=2.5", "2.6-3.0", "3.1-3.5", "3.6+"))
#'
#' irisNew$Petal.Length.Cat <- factor(cut(irisNew$Petal.Length, c(0, 2, 4, 6, Inf)),
#'                                  labels = c("<=2.0", "2.1-4.0", "4.1-6.0", "6.0+"))
#'
#' irisNew$Petal.Width.Cat <- factor(cut(irisNew$Petal.Width, c(0, 1, 2, Inf)),
#'                                labels = c("<=1.0", "1.1-2.0", "2.1+"))
#' 
#' ## Using NB to predict if the species is virginica
#' ## (training and predicting on same dataset)
#' pred <- performNB(irisNew, irisNew, obsIDVar = "id", goldStdVar = "spVirginica",
#' covariates = c("Sepal.Length.Cat", "Sepal.Width.Cat",
#'                "Petal.Length.Cat", "Petal.Width.Cat"), l = 1)
#' irisResults <- merge(irisNew, pred$probabilities, by = "id")
#' tapply(irisResults$p, irisResults$Species, summary)
#' 
#' @export



performNB <- function(training, prediction, obsIDVar, goldStdVar, 
                      covariates, l = 1){
  
  #### Checking variable names ####
  
  #Checking that the named variables are in the dataframe
  if(!obsIDVar %in% names(training) | !obsIDVar %in% names(prediction)){
    stop(paste0(obsIDVar, " is not in the dataframe."))
  }
  if(!goldStdVar %in% names(training)){
    stop(paste0(goldStdVar, " is not in the dataframe."))
  }
  #Checking that the covariates are in the dataframe
  covarTestT <- c(covariates %in% names(training), covariates %in% names(prediction))
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
  covarDf <- prediction[, covariates]
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
    probs <- prediction
    probs$p <- NA
    if(!"p" %in% names(training)){training$p <- NA}
    probs <- dplyr::bind_rows(probs[, c("p", obsIDVar)], training[, c("p", obsIDVar)])
    coeff <- NULL
    
    warning("No events or no non-events in training set")
    
  }else{
    
    #Setting up weight table
    #(currently irrelevant but would allow for weighting in subsequent versions)
    varTable <- as.data.frame(covariates)
    varTable$variable <- as.character(covariates)
    varTable$weight <- 1
    
    #Creating the results dataframe which is a copy of the prediction dataframe
    results <- prediction
    #Finding proportion of events/non-events in the trianing dataset
    classTab <- prop.table(table(training[, goldStdVar]) + l)
    
    #Initializing the coefficient dataframe
    coeff <- data.frame("level" = character(), "odds" = numeric())
    
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
      
      #Saving the odds: P(Outcome = TRUE|X=x, W) / P(Outcome = FALSE|X=x, W)
      odds <- Tab[, 2] / Tab[, 1]
      or <- odds/odds[1]
      num <- W * table(training[, Var], training[, goldStdVar]) + l
      se <- NA
      for(i in 1:(nrow(num) - 1)){
        numTab <- num[c(1, i+1), ]
        se <- c(se, sqrt(sum(1/numTab)))
      }
      level <- paste(Var, names(or), sep = ":")
      cTemp <- cbind.data.frame(level, or, se, stringsAsFactors = FALSE)

      coeff <- dplyr::bind_rows(coeff, cTemp)
    }
    
    #Calculating numerator and denominator for the probability calculation
    results$event <- apply(results[, grepl("_T", names(results))], 1, prod) * classTab[2] 
    results$nonevent <- apply(results[, grepl("_F", names(results))], 1, prod) * classTab[1]
    
    #Calculating probability of event
    probs <- results
    probs$p <- probs$event / (probs$event + probs$nonevent)
    
    #Combining training and prediction datasets
    if(!"p" %in% names(training)){training$p <- NA}
    probs <- dplyr::bind_rows(probs[, c(obsIDVar, "p")], training[, c(obsIDVar, "p")])
  }
  
  return(list("probabilities" = probs, "estimates" = coeff))
}

