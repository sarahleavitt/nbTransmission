# nbTransmission

[![Build Status](https://travis-ci.com/sarahleavitt/nbTransmission.png)](https://travis-ci.com/sarahleavitt/nbTransmission) [![codecov](https://codecov.io/gh/sarahleavitt/nbTransmission/branch/master/graph/badge.svg)](https://codecov.io/gh/sarahleavitt/nbTransmission)

For documentation and a tutorial see: https://sarahleavitt.github.io/nbTransmission/

## Introduction

This package is a group of functions used in infectious diseases analysis.
It implements an algorithm to calculate relative transmission probabilities between
cases in an infectious disease outbreak or cluster using naive Bayes. It also
contians various functions to use these probabilities to estimate
transmission paramaters such as the serial interval and reproductive number as well as
estimate the contribution of covariates to the probabilities and visualize the results.  

The ideal use of this package is for infectious disease dataset with metadata on the
majority of cases but more informative data such as contact tracing or pathogen whole
genome sequencing (WGS) on only a subset of cases. The package's algorithm allows
a researcher to infer transmission patterns among all cases and not just those
with the WGS or contact investigation data.  

Naive Bayes is a simple machine learning method that uses Bayes rule to estimate 
the probability of an outcome in a prediction dataset given a set of covariates 
from the observed frequencies in a training dataset. In this application, the outcome
is whether a pair is linked by direct transmission and the covariates could be spatial,
clinical, demographic, and temporal characteristics of the pairs of cases. A subset 
of cases with pathogen WGS or contact investigation data are used to create a
training dataset of probable transmission links and non/links and the relative probability
of a transmission link is estimated for all pairs.

For a more formal discussion of the theory behind and usage of this method, see the following paper:

Sarah V Leavitt, Robyn S Lee, Paola Sebastiani, C Robert Horsburgh, Jr,, Helen E Jenkins, Laura F White, Estimating the relative probability of direct transmission between infectious disease patients, International Journal of Epidemiology, dyaa031, https://doi.org/10.1093/ije/dyaa031
 
## Installation

You can install nbTransmission in R using the following command:

`devtools::install_github('https://github.com/sarahleavitt/nbTransmission.git')`


## Tutorial
Included in the package is a vingette that walks through how to use this method to analyze an infectious disease outbreak using the simulated datasets also included in this package. You can also access the tutorial at: https://sarahleavitt.github.io/nbTransmission/articles/nbTransmission-vignette.html

You could also install the vingette with the package using the following command (this may take a bit longer):

`devtools::install_github('https://github.com/sarahleavitt/nbTransmission.git', build_vignettes = TRUE)`

If you need assistance using nbTransmission, you can email sv1205@bu.edu.
