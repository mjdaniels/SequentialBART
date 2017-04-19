
<!-- README.md is generated from README.Rmd. Please edit that file -->
Description of the package
==========================

This program is used for imputing missing covariates by the 'sequential BART' approach. Sequential BART is a flexible Bayesian nonparametric approach to impute the missing covariates which involves factoring the joint distribution of the covariates with missingness into a set of sequential conditionals and applying Bayesian additive regression trees (BART) to model each of these univariate conditionals. Package provides a function, `seqBART()`, which computes and returns the imputed values.

Installation
============

`devtools::install_github("mjdaniels/SequentialBART")`

Main components of the package
==============================

The pacakge provides a function, `seqBART()`, to run the sequential BART model to find the missing covariates. The function takes as arguments

1.  x, Covariates having the missing values.

2.  y, Response Variable.

3.  x.type, a vector indicating the type of covariates (0=binary, 1=continuous)

4.  y.ype, the type of response and the inference regression model used for imputation, representing the type of missingness of the covaraites. It can take 5 values: y.type=0 for no response, y.type=1 for continuous response using linear regression for imputation, y.type=2 for binary response using logistic regression for imputation, y.type=3 for continuous response using BART for imputation, y.type=4 for binary response using BART probit for imputation.

Rest of the arguments are standard arguments for BART; Descriptions and defaults are provided.

Example
=======

`sbart::seqBART(x=Xcovariates, y=Response, x.type=datatypeValues, y.type=1)`
