#Sys.setenv("R_TESTS" = "")

library(testthat)
library(bartpkg1)

test_check("bartpkg1")
