#' Trait-aligned projection in PCoA space (multiple algorithms)
#'
#' For a given \eqn{n \times n} beta-diversity matrix \eqn{D} and a continuous
#' trait vector \eqn{y}, compute a **one-dimensional score**
#' \eqn{\hat y = E\hat\beta} from the first \code{ndim} PCoA axes \eqn{E}.
#'  
#' The score is used *only* for colouring the PCoA plot; the user can choose
#' one of three projection rules:
#'
#' \itemize{
#'   \item \strong{none}   — skip projection entirely and colour points by the
#'         raw trait.
#'   \item \strong{orthogonal}   — ordinary least-squares direction
#'         \eqn{\hat\beta = (E^\top E)^{-1}E^\top y}; returns adj.\;R\eqn{^2}
#'         and global F-test \eqn{p}.
#'   \item \strong{ordered}   — solve
#'         \deqn{\min_\beta \|y - E\beta\|_2^2 \quad
#'               \text{s.t. } (E_i-E_j)\beta > 0 \text{ if } y_i>y_j,\;
#'               (E_i-E_j)\beta = 0 \text{ if } y_i=y_j,}
#'         using \pkg{CVXR}.  
#'         Constraints are applied only to \emph{consecutive} trait-ordered
#'         samples for efficiency.  A non-negative \eqn{\ell_1} penalty
#'         (\code{lambda}) can be added to stabilise the solution when
#'         \eqn{ndim \gg n}.
#' }
#'
#' @param distance_matrix  A symmetric distance matrix (class \code{dist} or
#'   numeric square matrix) of beta diversities.
#' @param trait            Numeric vector of length \eqn{n} with the continuous
#'   phenotype.
#' @param ndim             Number of PCoA axes retained (default \code{5}).
#' @param figure           Logical; if \code{TRUE} (default) draw a PCoA1–PCoA2
#'   scatter coloured by the chosen projection or by the trait.
#' @param projection       Character string: \code{"none"}, \code{"orthogonal"}
#'   or \code{"ordered"} (see details above).
#' @param lambda           Optional non-negative LASSO weight used \emph{only}
#'   when \code{projection = "ordered"}; default \code{0}.
#'
#' @return A list with the same structure as the original implementation:
#'   \describe{
#'     \item{\code{projection_scores}}{Length-\eqn{n} vector used for colouring
#'           the plot (raw trait if \code{projection = "none"}).}
#'     \item{\code{gradient_direction}}{\eqn{\hat\beta}; \code{NA} when
#'           \code{projection = "none"}.}
#'     \item{\code{PCoA_coordinates}}{\eqn{n \times ndim} matrix of PCoA axes.}
#'     \item{\code{adjusted_R2}, \code{p_value}}{Model diagnostics returned only
#'           for the orthogonal fit; otherwise \code{NA}.}
#'     \item{\code{plot}}{A \pkg{ggplot2} object (or \code{NULL} if
#'           \code{figure = FALSE}).}
#'     \item{\code{variance_explained}}{Vector of eigenvalue proportions.}
#'   }
#' @importFrom stats cmdscale lm pf
#' @importFrom ggplot2 ggplot aes geom_point labs theme_classic
#' @importFrom viridis scale_color_viridis
#' @importFrom CVXR Variable sum_squares p_norm Problem Minimize solve
#' @export
gradient_projection <- function(distance_matrix,
                                trait,
                                ndim       = 5,
                                figure     = TRUE,
                                projection = c("orthogonal", "none", "ordered"),
                                lambda     = 0) {
  
  projection <- match.arg(projection)
  
  ## ---- dependencies ------------------------------------------------------
  if (!requireNamespace("vegan", quietly = TRUE))
    stop("Please install the 'vegan' package.")
  if (projection == "ordered" && !requireNamespace("CVXR", quietly = TRUE))
    stop("Projection 'ordered' requires the 'CVXR' package.")
  if (figure && !requireNamespace("ggplot2", quietly = TRUE))
    stop("Please install the 'ggplot2' package.")
  if (figure && !requireNamespace("viridis", quietly = TRUE))
    stop("Please install the 'viridis' package.")
  
  ## ---- input checks ------------------------------------------------------
  if (!is.matrix(distance_matrix) && !inherits(distance_matrix, "dist"))
    stop("distance_matrix must be matrix or 'dist'.")
  if (inherits(distance_matrix, "dist"))
    distance_matrix <- as.matrix(distance_matrix)
  
  n <- length(trait)
  if (nrow(distance_matrix) != n)
    stop("trait length must match number of samples.")
  
  ## ---- PCoA --------------------------------------------------------------
  pcoa_obj <- cmdscale(distance_matrix, k = ndim, eig = TRUE)
  E        <- pcoa_obj$points
  eigvals  <- pcoa_obj$eig
  var_expl <- eigvals / sum(eigvals[eigvals > 0])
  
  ## ---- choose projection algorithm --------------------------------------
  v_hat <- rep(NA_real_, ndim)
  proj  <- trait                      # default if "none"
  adj_r2 <- p_val <- NA_real_
  
  if (projection == "orthogonal") {
    v_hat <- solve(t(E) %*% E, t(E) %*% trait)
    proj  <- as.vector(E %*% v_hat)
    
    lm_fit <- lm(trait ~ E)
    fstat  <- summary(lm_fit)$fstatistic
    p_val  <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
    adj_r2 <- summary(lm_fit)$adj.r.squared
    
  } else if (projection == "ordered") {
    
    ## build monotone constraints only between consecutive sorted samples
    ord   <- order(trait)
    E_ord <- E[ord, , drop = FALSE]
    
    beta  <- CVXR::Variable(ncol(E))
    obj   <- CVXR::sum_squares(trait - E %*% beta)
    if (lambda > 0) obj <- obj + lambda * CVXR::p_norm(beta, 1)
    
    con <- list()
    for (i in 1:(n - 1)) {
      diffE <- E_ord[i + 1, ] - E_ord[i, ]
      if (trait[ord[i + 1]] > trait[ord[i]])
        con[[length(con) + 1]] <- t(diffE) %*% beta >= 1e-6
      else
        con[[length(con) + 1]] <- t(diffE) %*% beta == 0
    }
    
    prob <- CVXR::Problem(CVXR::Minimize(obj), con)
    sol  <- CVXR::solve(prob)
    if (sol$status != "optimal")
      warning("CVXR solver status: ", sol$status)
    
    v_hat <- as.numeric(sol$getValue(beta))
    proj  <- as.vector(E %*% v_hat)
    ## R² not meaningful here => keep NA
  }
  
  ## ---- plotting ----------------------------------------------------------
  p <- NULL
  if (figure) {
    df <- data.frame(PCoA1 = E[, 1], PCoA2 = E[, 2], Color = proj)
    
    title_txt <- switch(projection,
                        none       = "PCoA coloured by trait",
                        orthogonal = "PCoA coloured by OLS projection",
                        ordered    = "PCoA coloured by monotone projection")
    
    p <- ggplot2::ggplot(df,
                         ggplot2::aes(PCoA1, PCoA2, color = Color)) +
      ggplot2::geom_point(size = 3, alpha = 0.9) +
      viridis::scale_color_viridis(option = "plasma", direction = -1,
                                   name = "Colour scale") +
      ggplot2::labs(
        title = title_txt,
        x     = sprintf("PCoA1 (%.1f%%)", var_expl[1] * 100),
        y     = sprintf("PCoA2 (%.1f%%)", var_expl[2] * 100)
      ) +
      ggplot2::theme_classic(base_size = 14)
    print(p)
  }
  
  ## ---- return ------------------------------------------------------------
  list(
    projection_scores  = proj,
    gradient_direction = v_hat,
    PCoA_coordinates   = E,
    adjusted_R2        = adj_r2,
    p_value            = p_val,
    method             = projection,
    plot               = p,
    variance_explained = var_expl
  )
}