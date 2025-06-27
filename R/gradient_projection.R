#' Trait-Aligned Gradient Projection in PCoA Space
#'
#' Compute the 1-D projection of a continuous trait onto a host-community
#' PCoA embedding (beta-diversity space), following the formulation
#' \eqn{\hat v = (E^\top E)^{-1} E^\top y}.  The function returns
#' sample-level projection scores, the gradient direction vector,
#' adjusted \eqn{R^{2}}, and (optionally) a publication-quality plot.
#'
#' @param distance_matrix Either a numeric \code{dist} object or an
#'   \eqn{n \times n} symmetric distance matrix of community dissimilarities.
#' @param trait Numeric vector of length \eqn{n} giving the continuous
#'   phenotype/trait.
#' @param ndim Integer. Number of PCoA axes to retain (default \code{5}).
#' @param figure Logical. If \code{TRUE} (default) a ggplot2/viridis figure is
#'   printed and returned.
#'
#' @return A list with elements:
#' \describe{
#'   \item{projection_scores}{Numeric vector of length \eqn{n}; the trait-aligned
#'         scores \eqn{\hat y_i = E_i^\top \hat v}.}
#'   \item{gradient_direction}{Numeric vector of length \code{ndim}; the
#'         direction \eqn{\hat v}.}
#'   \item{PCoA_coordinates}{\eqn{n \times ndim} matrix of PCoA scores.}
#'   \item{adjusted_R2}{Adjusted \eqn{R^{2}} from the regression
#'         \code{trait ~ E}.}
#'   \item{p_value}{Global F-test p-value testing trait–embedding association.}
#'   \item{plot}{\code{ggplot} object (or \code{NULL} if \code{figure = FALSE}).}
#'   \item{variance_explained}{Vector of eigenvalue proportions.}
#' }
#'
#' @examples
#' \dontrun{
#'   library(vegan)
#'   set.seed(1)
#'   otu   <- matrix(abs(rnorm(100 * 30)), nrow = 100)
#'   trait <- rnorm(100)
#'   dist_mat <- vegdist(otu, method = "bray")
#'
#'   res <- gradient_projection(dist_mat, trait)
#'   res$plot            # interactive figure
#'   res$adjusted_R2
#' }
#'
#' @importFrom stats cmdscale dist lm pf
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal
#' @importFrom viridis scale_color_viridis
#' @export
gradient_projection <- function(distance_matrix, trait, ndim = 5, figure = TRUE) {
  # Load required packages
  if (!requireNamespace("vegan", quietly = TRUE)) stop("Please install the 'vegan' package.")
  if (figure && !requireNamespace("ggplot2", quietly = TRUE)) stop("Please install the 'ggplot2' package.")
  if (figure && !requireNamespace("viridis", quietly = TRUE)) stop("Please install 'viridis')")
  # Input checks
  if (!is.matrix(distance_matrix) && !inherits(distance_matrix, "dist")) {
    stop("distance_matrix must be a matrix or 'dist' object.")
  }
  
  n <- length(trait)
  if (inherits(distance_matrix, "dist")) {
    if (attr(distance_matrix, "Size") != n) stop("trait length must match number of samples.")
  } else {
    if (nrow(distance_matrix) != n) stop("trait length must match number of samples.")
  }
  
  # Step 1: PCoA with eigenvalues
  pcoa_obj <- cmdscale(distance_matrix, k = ndim, eig = TRUE)
  E <- pcoa_obj$points  # coordinates
  eigvals <- pcoa_obj$eig
  
  # Variance explained by first two axes
  var_explained <- eigvals / sum(eigvals[eigvals > 0])  # only positive eigenvalues
  var1 <- round(var_explained[1] * 100, 1)
  var2 <- round(var_explained[2] * 100, 1)
  
  # Step 2: Gradient direction
  v_hat <- solve(t(E) %*% E) %*% t(E) %*% trait
  y_hat <- as.vector(E %*% v_hat)
  
  # Step 3: Linear model fit and p-value
  lm_fit <- lm(trait ~ E)
  fstat <- summary(lm_fit)$fstatistic
  p_value <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  adj_r_squared <- summary(lm_fit)$adj.r.squared
  
  # Step 4: Optional figure
  p <- NULL
  if (figure) {
    library(ggplot2)
    library(viridis)
    
    plot_df <- data.frame(
      PCoA1 = E[, 1],
      PCoA2 = E[, 2],
      Trait = trait,
      Projection = y_hat
    )
    
    color_by <- if (p_value < 0.05) "Projection" else "Trait"
    plot_df$Color <- if (color_by == "Projection") y_hat else trait
    
    p <- ggplot2::ggplot(plot_df, aes(x = PCoA1, y = PCoA2, color = Color)) +
      ggplot2::geom_point(size = 3, alpha = 0.85) +
      viridis::scale_color_viridis(option = "plasma", direction = -1, name = color_by) +
      ggplot2::labs(
        title = "Gradient-Aligned Projection in PCoA Space",
        subtitle = if (color_by == "Projection") "Colored by projection scores" else "Colored by trait values",
        x = paste0("PCoA1 (", var1, "%)"),
        y = paste0("PCoA2 (", var2, "%)")
      ) +
      ggplot2::theme_minimal(base_size = 16) +
      ggplot2::theme(
        plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(size = 14, margin = margin(b = 10)),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    print(p)
  }
  
  
  return(list(
    projection_scores = y_hat,
    gradient_direction = v_hat,
    PCoA_coordinates = E,
    adjusted_R2 = adj_r_squared,
    p_value = p_value,
    plot = p,
    variance_explained = var_explained
  ))
}
