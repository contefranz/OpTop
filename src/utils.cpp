#include <RcppArmadillo.h>
#include "utils.h"

arma::mat normalize_columns(const arma::mat& x)
{
    arma::rowvec norms = arma::sqrt(arma::sum(x % x, 0));
    return x.each_row() / norms;
}
