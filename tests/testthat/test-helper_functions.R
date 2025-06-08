initChol.spam_old = function(T, D = 1){

  # Random initialization
  QHt_Matrix = build_Q(obs_sigma_t2 = abs(rnorm(T)),
                       evol_sigma_t2 = abs(rnorm(T)),
                       D = D)

  # And return the Cholesky piece:
  chQht_Matrix0 = chol.spam(as.spam.dgCMatrix(as(QHt_Matrix, "dgCMatrix")))

  chQht_Matrix0
}

test_that("new initChol works same as old", {
  set.seed(1)
  old_res <- initChol.spam_old(20)
  set.seed(1)
  res <- initChol_spam(20)

  expect_equal(as.matrix(old_res), as.matrix(res))
})
