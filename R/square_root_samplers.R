# Internal Gaussian recovery helpers. The precision expression is evaluated
# first; the square-root expression is built only when recovery is needed.
.btf_try_precision <- function(draw, fallback) {
  result <- tryCatch(withCallingHandlers({
    result <- force(draw)
    if(any(!is.finite(result))) stop("Nonfinite precision draw")
    result
  }, warning = function(w) stop(conditionMessage(w), call. = FALSE)),
  error = function(e) e)
  if(!inherits(result, "error")) return(result)

  result <- force(fallback)
  if(any(!is.finite(result))) stop("Nonfinite square-root recovery draw")
  .dsp_state$n_sqrt <- .dsp_state$n_sqrt + 1L
  result
}

.btf_chol <- function(Q) {
  R <- Matrix::chol(Q)
  diagonal <- Matrix::diag(R)
  if(any(!is.finite(diagonal)) || any(diagonal <= 0))
    stop("Invalid precision Cholesky diagonal")
  R
}

# Identity rows for the first D raw heights, followed by Dth differences.
btf_difference_operator <- function(n, D) {
  if(length(D) != 1L || is.na(D) || !D %in% 0:2 || n <= D)
    stop("Square-root BTF requires D in 0:2 and a state length greater than D")
  if(D == 0) return(Matrix::Diagonal(n))
  rows <- seq.int(D + 1L, n)
  Matrix::sparseMatrix(
    i = c(seq_len(D), rep(rows, each = D + 1L)),
    j = c(seq_len(D), as.numeric(t(outer(rows, -D:0, "+")))),
    x = c(rep(1, D), rep(choose(D, 0:D)*(-1)^(D:0), length(rows))),
    dims = c(n, n))
}

.btf_variances <- function(x, n) {
  if(!is.numeric(x) || !length(x) %in% c(1L, n) ||
     any(!is.finite(x)) || any(x <= 0))
    stop("Square-root BTF variances must be finite and positive, with length 1 or the state/observation length")
  rep_len(as.numeric(x), n)
}

# Draw N((A'A)^-1 A'z, (A'A)^-1) without forming the precision A'A.
# A must have full column rank. For A*S*P = U*R, solve
# R*theta = U'z + N(0,I), then return S*P*theta. Orthogonal transformations
# are applied implicitly; neither a dense U nor a matrix inverse is formed.
sample_gaussian_sqrt <- function(A, z, column_scale = TRUE, qr_order = 0L) {
  A <- as(as(as(Matrix::Matrix(A, sparse = TRUE), "generalMatrix"),
             "dMatrix"), "CsparseMatrix")
  z <- as.numeric(z)
  n <- ncol(A)
  if(n < 1L || nrow(A) < n || length(z) != nrow(A) ||
     any(!is.finite(A@x)) || any(!is.finite(z)))
    stop("Invalid square-root Gaussian operator or response")
  if(!is.logical(column_scale) || length(column_scale) != 1L || is.na(column_scale))
    stop("column_scale must be TRUE or FALSE")
  if(length(qr_order) != 1L || !is.finite(qr_order) || !qr_order %in% 0:3)
    stop("qr_order must be an integer in 0:3")

  scale <- rep(1, n)
  if(column_scale) {
    for(j in seq_len(n)) {
      first <- A@p[j] + 1L
      last <- A@p[j + 1L]
      if(first > last) stop("Square-root Gaussian operator is rank deficient")
      values <- A@x[seq.int(first, last)]
      maximum <- max(abs(values))
      if(maximum == 0) stop("Square-root Gaussian operator is rank deficient")
      norm <- maximum*sqrt(sum((values/maximum)^2))
      if(!is.finite(norm) || norm <= 0)
        stop("Invalid square-root Gaussian column norm")
      scale[j] <- 1/norm
    }
    A <- A %*% Matrix::Diagonal(x = scale)
  }

  factor <- Matrix::qr(A, order = as.integer(qr_order))
  R <- as(Matrix::qr.R(factor, complete = FALSE, backPermute = FALSE), "dtCMatrix")
  diagonal <- Matrix::diag(R)
  if(any(!is.finite(diagonal)) || any(diagonal == 0))
    stop("Square-root Gaussian operator is numerically rank deficient")
  projected <- as.numeric(Matrix::qr.qty(factor, z))[seq_len(n)]
  permuted <- as.numeric(Matrix::solve(R, projected + stats::rnorm(n)))
  permutation <- if(length(factor@q)) as.integer(factor@q) + 1L else seq_len(n)
  unpermuted <- numeric(n)
  unpermuted[permutation] <- permuted
  result <- matrix(scale*unpermuted, ncol = 1L)
  if(any(!is.finite(result))) stop("Nonfinite square-root Gaussian draw")
  result
}

sampleBTF_sqrt <- function(y, obs_sigma_t2, evol_sigma_t2, D,
                           prior_mean = NULL, linear_term = NULL,
                           predictor = NULL) {
  n <- length(y)
  v <- .btf_variances(obs_sigma_t2, n)
  w <- .btf_variances(evol_sigma_t2, n)
  if(is.null(predictor)) predictor <- rep(1, n)
  H <- Matrix::Diagonal(x = as.numeric(predictor)/sqrt(v))
  L <- btf_difference_operator(n, D)
  A <- rbind(H, Matrix::Diagonal(x = 1/sqrt(w)) %*% L)
  response <- if(is.null(linear_term)) y/sqrt(v) else linear_term*sqrt(v)
  result <- sample_gaussian_sqrt(A, c(response, numeric(n)))
  # Preserve sampleBTF's existing D>0 term b = V^-1*y + Q*prior_mean.
  # Its interpretation differs from a prior-only shift; correcting it is a
  # separate model change. Internal nonzero prior_mean calls currently use D=0.
  if(!is.null(prior_mean)) result <- result + as.numeric(prior_mean)
  result
}

sampleBTF_reg_sqrt <- function(y, X, obs_sigma_t2, evol_sigma_t2, D) {
  n <- nrow(X)
  p <- ncol(X)
  v <- .btf_variances(obs_sigma_t2, n)
  w <- .btf_variances(as.numeric(t(evol_sigma_t2)), n*p)
  # Time-major ordering: beta[1,1:p], beta[2,1:p], ... .
  H <- Matrix::sparseMatrix(i = rep(seq_len(n), each = p),
    j = seq_len(n*p), x = as.numeric(t(X/sqrt(v))), dims = c(n, n*p))
  L <- kronecker(btf_difference_operator(n, D), Matrix::Diagonal(p))
  A <- rbind(H, Matrix::Diagonal(x = 1/sqrt(w)) %*% L)
  matrix(sample_gaussian_sqrt(A, c(y/sqrt(v), numeric(n*p))),
         nrow = n, byrow = TRUE)
}

sampleBTF_bspline_sqrt <- function(y, X, obs_sigma2, evol_sigma_t2, D) {
  v <- .btf_variances(obs_sigma2, 1L)
  A <- Matrix::Matrix(X/sqrt(v), sparse = TRUE)
  z <- y/sqrt(v)
  # Preserve the current D=0 precision, which contains X'X/v only.
  # Adding its omitted independent prior would be a separate correction.
  if(D != 0) {
    p <- ncol(X)
    w <- .btf_variances(as.numeric(evol_sigma_t2), p)
    A <- rbind(A, Matrix::Diagonal(x = 1/sqrt(w)) %*% btf_difference_operator(p, D))
    z <- c(z, numeric(p))
  }
  sample_gaussian_sqrt(A, z)
}
