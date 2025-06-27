#' G-C-O Structural Decomposition
#'
#' Classifies each sample as \strong{G-type} (gradient-driven),
#' \strong{C-type} (cluster-driven), or \strong{O-type} (outlier) by comparing
#' neighbourhood structure in the original PCoA space to that in the 1-D
#' trait-aligned projection.
#'
#' \enumerate{
#'   \item For every sample the \code{k} nearest neighbours are found in
#'         PCoA (\eqn{\mathcal N_E}) and in projection space
#'         (\eqn{\mathcal N_G}).
#'   \item The local \emph{consistency score}
#'         \eqn{s_i = |\mathcal N_E \cap \mathcal N_G|/k} is computed.
#'   \item A sample is labelled:
#'     \describe{
#'       \item{G-type}{if \eqn{s_i > \tau_H} and not an outlier;}
#'       \item{C-type}{if \eqn{s_i < \tau_L} and not an outlier;}
#'       \item{O-type}{if the mean distance to its \code{k} neighbours exceeds
#'         the 95th percentile in \emph{either} space.}
#'     }
#' }
#'
#' @param pcoa_mat Numeric \eqn{n \times d} matrix of PCoA coordinates.
#' @param projection_scores Numeric vector (length \eqn{n}) of 1-D projection scores.
#' @param k Integer, neighbourhood size (default \code{10}).
#' @param tau_H High threshold for gradient consistency (default \code{0.8}).
#' @param tau_L Low  threshold for cluster consistency (default \code{0.2}).
#'
#' @return A data frame with one row per sample:
#' \describe{
#'   \item{\code{consistency_score}}{\eqn{s_i}.}
#'   \item{\code{outlier}}{Logical flag.}
#'   \item{\code{class}}{Factor with levels “G-type”, “C-type”, “O-type”, “Intermediate”.}
#' }
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   coords <- matrix(rnorm(40 * 2), ncol = 2)
#'   proj   <- coords %*% c(1, 0.5) + rnorm(40)
#'   gco    <- gco_decomposition(coords, proj, k = 8)
#'   table(gco$class)
#' }
#'
#' @importFrom stats dist
#' @export
gco_decomposition <- function(pcoa_mat, projection_scores, k = 10, tau_H = 0.8, tau_L = 0.2) {
  n <- nrow(pcoa_mat)
  
  # Distance in PCoA space
  dist_E <- as.matrix(dist(pcoa_mat))
  # Distance in gradient projection space (1D)
  dist_G <- as.matrix(dist(matrix(projection_scores, ncol = 1)))
  
  # Neighbor sets (k nearest, exclude self)
  get_kNN <- function(D) {
    apply(D, 1, function(row) {
      order(row)[2:(k + 1)]
    })
  }
  
  kNN_E <- get_kNN(dist_E)  # indices
  kNN_G <- get_kNN(dist_G)
  
  # Compute consistency score
  s_vec <- numeric(n)
  for (i in 1:n) {
    s_vec[i] <- length(intersect(kNN_E[, i], kNN_G[, i])) / k
  }
  
  # Define outliers based on average kNN distance > 95th percentile in either space
  kNN_avg_dist_E <- sapply(1:n, function(i) mean(dist_E[i, kNN_E[, i]]))
  kNN_avg_dist_G <- sapply(1:n, function(i) mean(dist_G[i, kNN_G[, i]]))
  
  threshold_E <- quantile(kNN_avg_dist_E, 0.95)
  threshold_G <- quantile(kNN_avg_dist_G, 0.95)
  
  outlier_flag <- (kNN_avg_dist_E > threshold_E) | (kNN_avg_dist_G > threshold_G)
  
  # Final classification
  class_label <- rep("Intermediate", n)
  class_label[s_vec > tau_H & !outlier_flag] <- "G-type"
  class_label[s_vec < tau_L & !outlier_flag] <- "C-type"
  class_label[outlier_flag] <- "O-type"
  
  return(data.frame(
    consistency_score = s_vec,
    outlier = outlier_flag,
    class = class_label
  ))
}
