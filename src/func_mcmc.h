// //---------------------------------
// //     ***** RCPP i.e. include<Rcpp.h> was Dropped by me
// //--------------------------------


// #include <Rcpp.h>
// using namespace Rcpp;

#include "rng.h"
#include "clasinit.h"
#include "funs.h"
#include "bd.h"

// //---------------------------------
// // function to get each mcmc run
// //--------------------------------
//
// //' run the mcmc functionality on data
// //' @param init object
// //' @export
// //'
// // [[Rcpp::export]]

void mcmc(init& ip_initial, RNG gen,int* vartype, size_t m, size_t n, double lambda, double nu);
