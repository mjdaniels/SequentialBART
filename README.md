
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

bartpkg1::seqBARTfunc(xx=xExample, yy=y1Example, datatype=datatypeExample, type=1)

``` r
# Commenting for now!!

# bartpkg1::serBARTfunc(xx=xExample, yy=y2Example, datatype=datatypeExample, type=2)

# only if the xx is present as a .rda and the pacakge is build once, then the above command works.

# It wont work if the xx was just in the sysdata.rda, while test and check were fine with sysdata.rda only! DONT KNOW WHY> SOME THING TO DO WITH LAZY LOADING as the build shows extra line 'moving datasets to lazyload DB' in the working case > BETTER TO HAVE THE README generted and remove those individual.rda files later so that you dont have to make them available to the user in the data() command.
```
