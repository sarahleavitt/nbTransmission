#' Performs naive bayes classification
#'
#' The function \code{performNB_BS} Calculates the posterior probabilities of a dichotomous class
#' variable given a set of covariates using Bayes rule and a bootstrapped adjusted odds ratio via logistic
#' regression.
#'
#' The main purpose of this function is to be used by \code{\link{nbProbabilities_BS}} to
#' estimate the relative transmission probability between individuals in an infectious
#' disease outbreak. However, it can be used more generally to estimate the probability
#' of any dichotomous outcome given a set of categorical covariates and adjusted odds ratios of such
#' dichotomous outcome.
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
#'  the outcome of interest.)
#' @param covariates A character vector containing the covariate variable names.
#' All covariates need to be categorical factor variables.
#' @param l Laplace smoothing parameter that is added to each cell
#' (a value of 0 indicates no smoothing).
#' @param nBS Number of bootstrap samples to run in each cross-validation fold/iteration (default is 100)
#' @param pSampled Proportion of unlinked cases to include in bootstrap sample (default is 1, i.e.a true bootstrap)
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
#'        \item \code{est} - the log odds ratio for this covariate and level
#'        \item \code{se} - the standard error of the log odds ratio
#'      }
#' }
#' @seealso \code{\link{nbProbabilities_BS}}
#'
#' @export

performNB_BS <- function(training, prediction, obsIDVar, goldStdVar,
                         covariates, l = 1,
                         nBS = 100, pSampled = 1){

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

  }

  else{

    #Setting up weight table
    #(currently irrelevant but would allow for weighting in subsequent versions)
    varTable <- as.data.frame(covariates)
    varTable$variable <- as.character(covariates)

    # ADD IN ATTRIBUTE WEIGHT
    varTable$weight <- 1

    #Creating the results dataframe which is a copy of the prediction dataframe
    results <- prediction

    #Finding proportion of events/non-events in the training dataset (prior -- P(c))
    classTab <- prop.table(table(training[, goldStdVar]) + l)

    #Initializing the coefficient dataframe
    coeff <- data.frame("level" = character(), "est" = numeric(),
                        "se" = numeric(), "Pr_F" = numeric(), "Pr_T" = numeric(), stringsAsFactors = FALSE)


    #Looping through all covariates and finding frequencies in training data
    for(i in 1:length(covariates)){

      #Extracting covariate name
      Var <- covariates[i]
      #Determining the weight for that covariate
      W <- varTable[varTable$variable == Var, "weight"]
      #Creating a table with proportions in each level of the covariate from training data
      Tab <- prop.table(W * table(training[, Var], training[, goldStdVar]) + l, 2)

      #Calculating P(X|Outcome = TRUE, W) and P(X|Outcome = FALSE, W) (P(a_j|c))
      #First creating a variable with the value of this covariate
      results$covariate <- results[, Var]
      #1 if missing, otherwise it is the proportion in the second column and the row
      #associated with the factor level of the covariate value.
      results[, paste0(Var, "_T")] <- ifelse(is.na(results$covariate), 1,
                                             Tab[as.numeric(results$covariate), 2]) ^ W
      results[, paste0(Var, "_F")] <- ifelse(is.na(results$covariate), 1,
                                             Tab[as.numeric(results$covariate), 1]) ^ W
    }
    if(pSampled == 1){

      bs_out <- list()

      for (i in 1:nBS) {
        # nSampled <- nrow(trainingF)*pSampled
        bs_samp <- sample_n(training, nrow(training), replace = T)
        # bs_samp <- rbind(trainingF_samp, trainingT)

        suppressWarnings(mylogit <- glm(reformulate(covariates, goldStdVar), data = bs_samp, family = "binomial"))
        tidy_fit <- tidy(mylogit)
        tidy_out <- data.frame(tidy_fit$term, tidy_fit$estimate)


        bs_out[[i]] <- tidy_out
      }
      # print("Running full bootstrap")
    }
    else{
      trainingT <- training %>% filter(linked == T)
      trainingF <- training %>% filter(linked == F)

      bs_out <- list()

      for (i in 1:nBS) {
        nSampled <- nrow(trainingF)*pSampled
        trainingF_samp <- sample_n(trainingF, nSampled, replace = T)
        bs_samp <- rbind(trainingF_samp, trainingT)

        suppressWarnings(mylogit <- glm(reformulate(covariates, goldStdVar), data = bs_samp, family = "binomial"))
        tidy_fit <- tidy(mylogit)
        tidy_out <- data.frame(tidy_fit$term, tidy_fit$estimate)


        bs_out[[i]] <- tidy_out


      }
      # print("Running smaller sample bootstrap")
    }

    bs_out1 <- bind_rows(bs_out)

    bs_out_df <- bs_out1 %>%
      group_by(tidy_fit.term) %>%
      summarise(bs_est = mean(tidy_fit.estimate),
                bs_sd = sd(tidy_fit.estimate)) %>%
      rename(level = tidy_fit.term)

    #Calculating numerator and denominator for the probability calculation
    results$event <- apply(results[, grepl("_T", names(results))], 1, prod) * classTab[2]
    results$nonevent <- apply(results[, grepl("_F", names(results))], 1, prod) * classTab[1]

    #Calculating probability of event
    probs0 <- results
    probs0$p <- probs0$event / (probs0$event + probs0$nonevent)

    #Combining training and prediction datasets
    if(!"p" %in% names(training)){training$p <- NA}
    probs <- dplyr::bind_rows(probs0[, c(obsIDVar, "p")], training[, c(obsIDVar, "p")])
  }

  return(list("probabilities" = probs, "estimates" = bs_out_df))
}
