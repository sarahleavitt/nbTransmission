
#' Pair-level Outbreak Dataset
#' 
#' A dataset of all pairs (unordered) of the outbreak of 100 individuals 
#' described in \code{\link{indData}} with SNP distance between each pair. Genomes were
#' simulated using the package phangorn from the phylogentic tree created
#' during the outbreak simulation. 
#' 
#' @format A data frame with 9900 rows and 17 variables:
#' \describe{
#' \item{pairID}{A pair-level ID variable (the individual IDs separated by an '_').}
#' \item{individualID.1}{The ID of the potential "infector".}
#' \item{individualID.2}{The ID of the potential "infectee".}
#' \item{transmission}{Did individual.1 infect individual.2.}
#' \item{snpDist}{What is the number of SNPs between the individuals.}
#' \item{infectionDate.1}{The date and time of infection of individualID.1.}
#' \item{infectionDate.2}{The date and time of infection of individualID.2.}
#' \item{sampleDate.1}{The date and time of sampling of individualID.1.}
#' \item{sampleDate.2}{The date and time of sampling of individualID.2.}
#' \item{sampleDiff}{The number of days between sampleDate.1 and sampleDate.2.}
#' \item{infectionDiff}{The number of days between infectionDate.1 and infectionDate.2.}
#' \item{infectionDiffY}{The number of years between infectionDate.1 and infectionDate.2.}
#' \item{timeCat}{A categorical representation of infectionDiff: <1y, 1-2y, 2-3y, 3-4y, 4-5y, >5y.}
#' \item{Z1}{Pair-level covariate derived from X1: 1 if match, 0 if not match.}
#' \item{Z2}{Pair-level covariate derived from X2: 1 if match, 0 if not match.}
#' \item{Z3}{Pair-level covariate derived from X3: 1 if a-a, 2 if b-b, 3 if a-b, 4 if b-a.}
#' \item{Z4}{Pair-level covariate derived from X4: 1 if match, 2 if adjacent, 2 otherwise.}
#' }
"pairData"
