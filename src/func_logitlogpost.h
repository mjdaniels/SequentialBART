#include "Eigen/Dense"
#include "Eigen/Cholesky"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using namespace Eigen;

double logit_logpost(MatrixXd Y, MatrixXd X, MatrixXd beta);
