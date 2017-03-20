
<!-- README.md is generated from README.Rmd. Please edit that file -->
Description of the package
==========================

This program is used for imputing missing covariates by the 'sequential BART' approach. The codes are revised based on the 'parallel BART' code f (though parallel computation of BART was not used).

Installation
============

devtools::install\_github("mjdaniels/SequentialBART")

Main components of the package
==============================

The pacakge provides a function 'SeqBART' to run the sequential BART model to find the missing covariates. The function takes as arguments 1. X, Covariates having the missing values.

1.  Y, Response Variable.

2.  Datatype, representing the type of covariate, continuous or binary,

3.  Type, representing the type of missingness of the covaraites. It can take 3 values: 0 to represent covariates are MAR with MDM not depending on the response, and 1 or 2 to represent covariates are MAR with MDM depending on the response. If the response is continuous, use type=1 ( linear regression used for imputation), else if it is binary, use type=2 (logistic regression used for imputation).

Rest of the arguments are standard values for proper imputation. Defaults are provided.

Example
=======

bartpkg1::serBARTfunc(xx=xExample, yy=y1Example, datatype=datatypeExample, type=1)
