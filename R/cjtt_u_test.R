#' CJTT-U Trend Test for Beta-Dispersion vs. Continuous Trait
#'
#' Implements the continuous-trait version of the U-statistic trend test
#' proposed on slides 9–12 (“CJTT-U”).  It tests whether sample-level beta
#' dispersion
#' \eqn{d_i = \frac{1}{n-1}\sum_{j\neq i} D_{ij}}
#' shows a monotone trend with a numeric trait
#' by computing
#' \deqn{U = \frac{1}{\binom{n}{2}} \sum_{i<j}
#'        \operatorname{sign}(y_j - y_i)\,
#'        \operatorname{sign}(d_j - d_i).}
#'
#' Variance can be estimated analytically, via leave-one-out jackknife, or by
#' non-parametric bootstrap.  A permutation option is also provided.
#'
#' If a matrix of confounders is supplied, double residualisation à la
#' Frisch–Waugh–Lovell (FWL) is applied to both the trait and the distance
#' matrix using \eqn{P = I - Z\,(Z^\top Z)^{-1}Z^\top}.  
#' **Note:** residualising a distance matrix directly is an approximation;
#' for exact FWL you should residualise ordination **coordinates** first and
#' then recompute distances.
#'
#' @param dist_matrix Either a numeric \code{dist} object or an
#'   \eqn{n\times n} distance matrix.
#' @param trait Numeric vector of length \eqn{n}.
#' @param confounder Optional \eqn{n\times p} matrix of confounders
#'   (default \code{NULL}).
#' @param var_method Character.  Variance estimator:
#'   \code{"analytic"} (default), \code{"jackknife"}, or \code{"bootstrap"}.
#' @param permutation Logical.  If \code{TRUE}, a permutation p-value is
#'   returned instead of a z-test.
#' @param B Integer.  Number of bootstrap or permutation replicates
#'   (default \code{1000}).
#' @param seed Integer for reproducibility.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{U_statistic}}{Value of U in \([-1,1]\).}
#'   \item{\code{variance}}{Estimated variance (NA if permutation).}
#'   \item{\code{z_score}}{Standardised statistic (NA if permutation).}
#'   \item{\code{p_value}}{Two-sided p-value.}
#'   \item{\code{method}}{Variance method actually used.}
#' }
#'
#' @examples
#' \dontrun{
#'   library(vegan)
#'   set.seed(1)
#'   otu   <- matrix(abs(rnorm(40 * 25)), nrow = 40)
#'   trait <- seq(1, 40) + rnorm(40, 0, 2)
#'   D     <- vegdist(otu, "bray")
#'
#'   # analytic variance
#'   cjtt_u_test(D, trait)
#'
#'   # jackknife with a confounder (e.g., age)
#'   age <- rnorm(40, 50, 8)
#'   cjtt_u_test(D, trait, confounder = age, var_method = "jackknife")
#' }
#'
#' @importFrom stats dist pf pnorm
#' @importFrom utils combn
#' @export
cjtt_u_test <- function(dist_matrix, trait, confounder = NULL,
                        var_method = c("analytic", "jackknife", "bootstrap"),
                        permutation = FALSE, B = 1000, seed = 123) {
  set.seed(seed)
  var_method <- match.arg(var_method)
  n <- length(trait)
  
  # --- Input checks ---
  if (!is.matrix(dist_matrix) && !inherits(dist_matrix, "dist")) {
    stop("dist_matrix must be a matrix or 'dist' object.")
  }
  if (inherits(dist_matrix, "dist")) {
    dist_matrix <- as.matrix(dist_matrix)
  }
  if (nrow(dist_matrix) != n) {
    stop("Trait length must match number of samples in distance matrix.")
  }
  
  # --- Optional: FWL residualization ---
  if (!is.null(confounder)) {
    Z <- as.matrix(confounder)
    H <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
    P <- diag(n) - H
    trait <- as.numeric(P %*% trait)
    dist_matrix <- P %*% dist_matrix %*% P
  }
  
  # --- Step 1: Compute sample-level dispersion (distance to others) ---
  dispersion <- rowSums(dist_matrix) / (n - 1)
  
  # --- Step 2: Compute U-statistic ---
  combs <- combn(n, 2)
  signs <- apply(combs, 2, function(idx) {
    i <- idx[1]; j <- idx[2]
    sign(trait[j] - trait[i]) * sign(dispersion[j] - dispersion[i])
  })
  U <- mean(signs)
  
  # --- Step 3: Permutation test (if selected) ---
  if (permutation) {
    perm_U <- numeric(B)
    for (b in 1:B) {
      perm_trait <- sample(trait)
      signs_b <- apply(combs, 2, function(idx) {
        i <- idx[1]; j <- idx[2]
        sign(perm_trait[j] - perm_trait[i]) * sign(dispersion[j] - dispersion[i])
      })
      perm_U[b] <- mean(signs_b)
    }
    p_val <- mean(abs(perm_U) >= abs(U))
    return(list(
      U_statistic = U,
      variance = NA,
      z_score = NA,
      p_value = p_val,
      method = "permutation"
    ))
  }
  
  # --- Step 4: Variance estimation ---
  if (var_method == "analytic") {
    var_hat <- (4 * (n - 2)) / (9 * n * (n - 1))
  }
  
  if (var_method == "jackknife") {
    U_i <- numeric(n)
    for (i in 1:n) {
      idx <- setdiff(1:n, i)
      trait_i <- trait[idx]
      disp_i <- dispersion[idx]
      combs_i <- combn(n - 1, 2)
      signs_i <- apply(combs_i, 2, function(k) {
        a <- k[1]; b <- k[2]
        sign(trait_i[b] - trait_i[a]) * sign(disp_i[b] - disp_i[a])
      })
      U_i[i] <- mean(signs_i)
    }
    U_bar <- mean(U_i)
    var_hat <- (n - 1) / n * sum((U_i - U_bar)^2)
  }
  
  if (var_method == "bootstrap") {
    U_b <- numeric(B)
    for (b in 1:B) {
      idx <- sample(1:n, n, replace = TRUE)
      trait_b <- trait[idx]
      disp_b <- dispersion[idx]
      combs_b <- combn(n, 2)
      signs_b <- apply(combs_b, 2, function(idx2) {
        i <- idx2[1]; j <- idx2[2]
        sign(trait_b[j] - trait_b[i]) * sign(disp_b[j] - disp_b[i])
      })
      U_b[b] <- mean(signs_b)
    }
    var_hat <- var(U_b)
  }
  
  # --- Step 5: z-score and p-value ---
  z <- U / sqrt(var_hat)
  p_val <- 2 * (1 - pnorm(abs(z)))
  
  return(list(
    U_statistic = U,
    variance = var_hat,
    z_score = z,
    p_value = p_val,
    method = var_method
  ))
}
