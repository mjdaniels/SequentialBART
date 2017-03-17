// change 1: function to calcualte logposterior

#include "func_logitlogpost.h"



double logit_logpost(MatrixXd Y, MatrixXd X, MatrixXd beta)
{
  // likelihood
  VectorXd eta = X * beta;
  // MatrixXd p = 1.0 / (1.0 + exp(-eta));
  double loglike = 0.0;

  for (unsigned int i = 0; i < Y.rows(); ++ i)
    loglike += -Y(i) * log(1.0+exp(-eta(i))) - (1 - Y(i)) * log(1.0+exp(eta(i)));


  return (loglike );
}
