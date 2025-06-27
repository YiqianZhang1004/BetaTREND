#' Evaluate Projection–PCoA Consistency
#'
#' Computes three standard neighbourhood‐preservation diagnostics that quantify
#' how faithfully a 1-D trait-aligned projection reproduces structure in a
#' higher-dimensional PCoA space:
#'
#' \itemize{
#'   \item \strong{Trustworthiness} — proportion of projection-space neighbours
#'         that are also close in the original PCoA space.
#'   \item \strong{Continuity} — proportion of PCoA-space neighbours that remain
#'         neighbours in the projection.
#'   \item \strong{Mantel correlation} — Pearson/Spearman correlation between
#'         the two distance matrices.
#' }
#'
#' Neighbour ranks are shifted by \eqn{-1} so the self-distance has rank 0 and
#' the nearest neighbour has rank 1 (matching the definition in van der Maaten &
#' Hinton 2008).
#'
#' @param pcoa_mat \eqn{n \times d} numeric matrix of PCoA (or any ordination)
#'   coordinates.
#' @param projection_scores Numeric vector of length \eqn{n} containing the
#'   1-D trait-aligned projection scores (e.g.\ output of
#'   \code{\link{gradient_projection}}).
#' @param k Integer. The neighbourhood size used for trustworthiness and
#'   continuity (default \code{10}).
#' @param method Character string passed to \code{vegan::\link[vegan]{mantel}}
#'   (“pearson” [default] or “spearman”).
#'
#' @return A list with
#' \describe{
#'   \item{\code{trustworthiness}}{Value in [0,1]; 1 = perfect.}
#'   \item{\code{continuity}}{Value in [0,1]; 1 = perfect.}
#'   \item{\code{mantel_statistic}}{Mantel \eqn{r}.}
#'   \item{\code{mantel_p}}{Two-sided Mantel p-value.}
#' }
#'
#' @references
#' van der Maaten, L.J.P. & Hinton, G.E. (2008) Visualizing Data using
#' t-SNE. \emph{Journal of Machine Learning Research}, \strong{9}, 2579–2605.
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   coords <- matrix(rnorm(30 * 3), ncol = 3)
#'   proj   <- coords %*% c(1, 0.5, -0.2) + rnorm(30, 0, 0.1)
#'
#'   metrics <- evaluate_consistency_metrics(coords, proj, k = 5)
#'   metrics$trustworthiness
#'   metrics$continuity
#'   metrics$mantel_statistic
#' }
#'
#' @importFrom stats dist
#' @importFrom vegan mantel
#' @export
evaluate_consistency_metrics <- function(pcoa_mat, projection_scores, k = 10, method = "pearson") {
  if (!requireNamespace("vegan", quietly = TRUE)) stop("Please install the 'vegan' package.")
  
  n <- nrow(pcoa_mat)
  
  # Distance matrices
  dist_E <- as.matrix(dist(pcoa_mat))                             # community (PCoA) space
  dist_G <- as.matrix(dist(matrix(projection_scores, ncol = 1)))  # gradient projection space
  
  # Rank matrices for Trustworthiness/Continuity
  rank_E <- apply(dist_E, 1, function(r) rank(r, ties.method = "average") - 1)
  rank_G <- apply(dist_G, 1, function(r) rank(r, ties.method = "average") - 1)
  
  # --- 1. Trustworthiness (trait-space neighbors in community space) ---
  trust_sum <- 0
  for (i in 1:n) {
    g_neighbors <- order(dist_G[i, ])[2:(k + 1)]  # exclude self
    for (j in g_neighbors) {
      rank_ij <- rank_E[i, j]
      if (rank_ij > k) trust_sum <- trust_sum + (rank_ij - k)
    }
  }
  denom <- n * k * (2 * n - 3 * k - 1)
  trustworthiness <- 1 - (2 / denom) * trust_sum
  
  # --- 2. Continuity (community-space neighbors in gradient space) ---
  continuity_sum <- 0
  for (i in 1:n) {
    e_neighbors <- order(dist_E[i, ])[2:(k + 1)]
    for (j in e_neighbors) {
      rank_ij <- rank_G[i, j]
      if (rank_ij > k) continuity_sum <- continuity_sum + (rank_ij - k)
    }
  }
  continuity <- 1 - (2 / denom) * continuity_sum
  
  # --- 3. Mantel correlation (between dist_E and |ŷ_i - ŷ_j|) ---
  mantel_result <- vegan::mantel(dist_E, dist_G, method = method)
  
  return(list(
    trustworthiness = trustworthiness,
    continuity = continuity,
    mantel_statistic = mantel_result$statistic,
    mantel_p = mantel_result$signif
  ))
}
