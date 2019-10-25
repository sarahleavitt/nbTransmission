
#' Individual-level Outbreak Dataset
#' 
#' A dataset of an outbreak of 100 individuals starting from one case.
#' It was simulated with using TransPhylo with parameters resembling TB,
#'  starting from one case:
#' \itemize{
#' \item average reproductive number = 1.2
#' \item serial interval = gamma(1.2, 2) years
#' \item outbreak duration = 14 years
#' \item mutation rate 0.5 snps/genome/year
#' }
#' 
#' @format A data frame with 100 rows and 8 variables:
#' \describe{
#' \item{individualID}{An individual-level id for each case}
#' \item{infector}{The individualID of the true infector}
#' \item{infectionDate}{The date and time of infection}
#' \item{sampleDate}{The date and time of sampling}
#' \item{X1}{Covariate with 2 values: a, b (e.g. sex)}
#' \item{X2}{Covariate with 4 values: a, b, c, d (e.g. nationality)}
#' \item{X3}{Covariate with 2 values: a, b (e.g. homelessness)}
#' \item{X4}{Covariate with 10 values: a-j (e.g. county of residence)}
#' }
"indData"

