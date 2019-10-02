# nbTransmission

## Introduction

This package is a group of functions used in infectious diseases analysis.
It implements an algorithm to calculate relative transmission probabilities between
cases in an infectious disease outbreak or cluster using naive Bayes. It also
contians various functions to use these probabilities to estimate
transmission paramaters such as the serial interval and reproductive number.  

The ideal use of this package is for infectious disease dataset with metadata on the
majority of cases but more informative data such as contact tracing or pathogen whole
genome sequencing (WGS) on only a subset of cases. The packages' algorithm allows
a researcher to infer transmission patterns amongst all cases and not just those
with the discrimatory information.  

Naive Bayes is a simple machine learning method that uses Bayes rule to estimate 
the probability of an outcome in a prediction dataset given a set of covariates 
from the observed frequencies in a training dataset. The covariates could be spatial,
clinical, demographic, and temporal characteristics of the cases. Then a subset 
of cases with pathogen WGS or contact investigation data are used to create a
 training dataset of probable links and non/links.  
 
## Installation

You can install nbTransmission in R using the following command:

`devtools::install_github('https://github.com/sarahleavitt/nbTransmission.git')`


## Documentation

IN PROGRESS
There is a tutorial for nbTransmission at ... which is also within the R package
as a vingette. The reference manual is located at ....  

If you need assistance using nbTransmission, you can email sv1205@bu.edu.
