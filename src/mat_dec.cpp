#include <vector>
#include <math.h>
#include <stdio.h>
#include <random>
#include <limits>
#include <Rcpp.h>
#include <RcppEigen.h>

typedef Eigen::SparseMatrix<double> SpMat;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::VectorXd sample_mat(std::vector<int> row_ind, std::vector<int> col_ind,
              std::vector<double> val, int n, int num_entries,
              Eigen::VectorXd linht, Eigen::VectorXd rd, int D)
{
  Eigen::SparseMatrix<double> mat(n, n);
  if (D == 1) {
    mat.reserve(Eigen::VectorXi::Constant(n,3));
  } else{
    mat.reserve(Eigen::VectorXi::Constant(n,5));
  }

  for (int j = 0; j < num_entries; j++) {
    mat.insert(row_ind[j], col_ind[j]) = val[j];
  }

  Eigen::SimplicialLLT <Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> chol;
  chol.compute(mat);

  // SimplicialLLT does not report failure through its result. If mat is not
  // numerically positive definite, matrixL() and matrixU() are left in an
  // unspecified state and the solve below returns meaningless values. These
  // are not reliably detectable by the caller: they may be NaN, but are also
  // sometimes finite values far outside the plausible range of the states, and
  // they vary between runs on identical input. Returning them silently
  // causes problems downstream, so check the status and signal failure as NaN.
  if (chol.info() != Eigen::Success) {
    return Eigen::VectorXd::Constant(n, std::numeric_limits<double>::quiet_NaN());
  }

  Eigen::VectorXd tp = chol.matrixU().solve(chol.matrixL().solve(linht) + rd);

  return tp;
}









