#include <vector>
#include <math.h>
#include <stdio.h>
#include <random>
#include <chrono>
#include <RcppArmadillo.h>


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::uvec draw_indicators_generic(std::vector<double> res, std::vector<int> r, int n,
                             std::vector<double> w, std::vector<double> m,
                             std::vector<double> s, int J)
{
  // res_i = \ep_i, \ep_i \sim N(m_{r_i}, v_{r_i}).
  std::vector<double> lpmf(J);
  std::vector<double> pmf(J);
  std::vector<double> cdf(J);
  // std::default_random_engine generator;
  // generator.seed(std::chrono::system_clock::now().time_since_epoch().count());

  std::vector<double> lcoef(J);
  for (int j = 0; j < J; j++)
    lcoef[j] = log(w[j]) - log(s[j]);

  for (int i = 0; i < n; i++) {

    // Construct CDF.
    cdf[0] = 0.0;
    for (int j = 0; j < J; j++) {
      lpmf[j] = lcoef[j] - 0.5 * pow((res[i] - m[j]) / s[j], 2);
      pmf[j] = exp(lpmf[j]);
      cdf[j] = pmf[j] + (j > 0 ? cdf[j-1] : 0.0);
    }

    // Draw from discrete.
    // std::uniform_real_distribution<double> distribution(0.0, cdf[J-1]);
    // double U = distribution(generator);
    double U = Rcpp::runif(1, 0.0, cdf[J-1])[0];
    int k = 0;
    while(cdf[k] < U)
      k++;
    r[i] = k+1; // k in C indexing, k+1 in R indexing.
  }

  return(conv_to<arma::uvec>::from(r));
}






