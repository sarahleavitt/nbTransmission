#' Performs naive bayes classification
#'
#' The function \code{performNB} calculates the posterior probabilities of a dichotomous class
#' variable given a set of covariates using Bayes rule and either a univariate (default, \code{orType = "univariate"}
#' odds ratio or a bootstrapped adjusted odds ratio via logistic regression.
#'
#' The main purpose of this function is to be used by \code{\link{nbProbabilities}} to
#' estimate the relative transmission probability between individuals in an infectious
#' disease outbreak. However, it can be used more generally to estimate the probability
#' of any dichotomous outcome given a set of categorical covariates and adjusted odds ratios of such
#' dichotomous outcome.
#'
#' This function also generates odds ratios describing the associations between covariates in the training data
#' and outcome defined in the gold standard variable (\code{goldStdVar}) argument. Unadjusted odds ratios are the default.
#' These odds ratios are produced using contingency table methods. Adjusted odds ratios are calculated via bootstrapped
#' logistic regression to produce non-parametric standard errors. The bootstrap is controlled by parameters \code{nBS},
#' the number of bootstrap samples to run, and \code{pSampled}, the proportion of unlinked cases to include in the bootstrap
#' sample. \code{pSampled} is recommended only for large datasets in which it is computationally unfeasible to run a full
#' bootstrap. Sensitivity analyses should be run to determine an adequate value for \code{pSampled}.
#'
#' The function needs a training dataset with the outcome variable (\code{goldStdVar})
#' which is \code{TRUE} for those who have the value of interest and \code{FALSE}
#' for those who do not. The probability of having the outcome
#' (\code{<goldStdVar> = TRUE}) is predicted in the prediction dataset.
#'
#' @import stats
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
#' @param orType Takes value \code{"univariate"} or \code{"adjusted"}. \code{"univariate"} produces contingency table
#' odds ratios and \code{"adjusted"} produces adjusted odds ratios from a bootstrapped multivariable logistic regression.
#' @param nBS Number of bootstrap samples to run in each cross-validation fold/iteration (default is 100). Only
#' relevant when \code{orType = "univariate"}.
#' @param pSampled Proportion of unlinked cases to include in bootstrap sample (default is 1, i.e.a true
#' bootstrap). Only relevant when \code{orType = "univariate"}.
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
#' @seealso \code{\link{nbProbabilities}}
#'
#' @export

performNB <- function(training, prediction, obsIDVar, goldStdVar,
                         covariates, l = 1,
                         orType = "univariate",
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
    varTable$weight <- 1

    #Creating the results dataframe which is a copy of the prediction dataframe
    results <- prediction
    #Finding proportion of events/non-events in the trianing dataset
    classTab <- prop.table(table(training[, goldStdVar]) + l)

    #Initializing the coefficient dataframe
    coeff <- data.frame("level" = character(), "est" = numeric(),
                        "se" = numeric(), stringsAsFactors = FALSE)


    #Looping through all covariates and finding frequencies in training data
    # print(orType)
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

      if(orType == "univariate"){
        # print("type = uni")
      #Saving the odds: P(Outcome = TRUE|X=x, W) / P(Outcome = FALSE|X=x, W)
      odds <- Tab[, 2] / Tab[, 1]
      est <- log(odds/odds[1])
      num <- W * table(training[, Var], training[, goldStdVar]) + l
      se <- NA

      for(i in 1:(nrow(num) - 1)){
        numTab <- num[c(1, i+1), ]
        se <- c(se, sqrt(sum(1/numTab)))
      }

      level <- paste(Var, names(est), sep = ":")
      cTemp <- cbind.data.frame(level, est, se, stringsAsFactors = FALSE)

      coeff <- dplyr::bind_rows(coeff, cTemp)}
    }

    if(orType == "adjusted"){
    if(pSampled == 1){
    # print("type = adj")
      bs_out <- list()

      for (i in 1:nBS) {

        bs_samp <- dplyr::sample_n(training, nrow(training), replace = T)

        suppressWarnings(mylogit <- stats::glm(stats::reformulate(covariates, goldStdVar), data = bs_samp, family = "binomial"))
        tidy_fit <- broom::tidy(mylogit)
        tidy_out <- data.frame(tidy_fit$term, tidy_fit$estimate)


        bs_out[[i]] <- tidy_out
      }
      # print("Running full bootstrap")
    }
    else{
      trainingT <- training[training$linked == T,]
      trainingF <- training[training$linked == F,]

      bs_out <- list()

      for (i in 1:nBS) {
        nSampled <- nrow(trainingF)*pSampled
        trainingF_samp <- dplyr::sample_n(trainingF, nSampled, replace = T)
        bs_samp <- rbind(trainingF_samp, trainingT)

        suppressWarnings(mylogit <- glm(reformulate(covariates, goldStdVar), data = bs_samp, family = "binomial"))
        tidy_fit <- broom::tidy(mylogit)
        tidy_out <- data.frame(tidy_fit$term, tidy_fit$estimate)


        bs_out[[i]] <- tidy_out
      }
      # print("Running smaller sample bootstrap")
    }

    bs_out1 <- dplyr::bind_rows(bs_out)

    coeff <- by(bs_out1,
                 INDICES = list(bs_out1$tidy_fit.term),
                 FUN = function(x){
                   data.frame("level" = unique(x$tidy_fit.term),
                              est = mean(x$tidy_fit.estimate),
                              se = stats::sd(x$tidy_fit.estimate))
                 })
    coeff <- do.call(dplyr::bind_rows, coeff)
    }

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

  return(list("probabilities" = probs, "estimates" = coeff))
}
