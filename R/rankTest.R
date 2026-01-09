#' @title Rank-score test for quantile regression
#' @description
#' Performs the rank-score test for the covariate of interest
#' \code{X}, at the quantiles defined in \code{tau}, using a fitted quantile
#' regression model. The test evaluates the null hypothesis that the coefficient of \code{X} is equal to zero
#' against a two-sided alternative, at each specified quantile level.
#' Testing equality to a non-zero value is not yet implemented.
#'
#' @usage rankTest(mod, X, tau = NULL, full = FALSE, h = NULL, alpha = 0.05,
#' eps = c(1e-04,1e-04), B = "identity", error.distr = NULL, error.par = NULL)
#'
#' @param mod An object of class \code{rqs} returned by
#'   \code{\link[quantreg]{rq}}, representing the fitted quantile regression models.
#' @param X A string indicating the covariate of interest.
#' @param tau A numeric vector of quantiles of interest used in \code{mod}.
#' If \code{NULL} (default), all quantiles from the \code{mod} object are considered.
#' @param full Logical. If \code{TRUE}, the function returns the test statistics
#' and corresponding \eqn{p}-values for all intersection hypotheses containing
#' \code{tau}. If \code{FALSE} (default), only the results for the single hypotheses
#' are returned.
#' @param h Character string specifying the bandwidth selection method for
#' sparsity estimation. One of \code{"Hall-Sheather"} or \code{"Bofinger"}.
#' The default, \code{NULL}, uses the Hall-Sheather (1988) rule.
#' @param alpha A numeric value used for bandwidth estimation.
#' Following Koenker (2005), it is typically set equal to the desired significance level.
#' @param eps Numeric vector of length 2 specifying the relative accuracies
#'   requested for approximating the distribution of the test statistic using
#'   Imhof's (1961) method, "Computing the Distribution of Quadratic Forms in
#'   Normal Variables". The first component controls the accuracy of the initial
#'   (more accurate) numerical integration (using \code{\link[stats]{integrate}}),
#'   while the second component controls the accuracy of a subsequent, less
#'   accurate integration step (using \code{\link[pracma]{quadgk}}).
#' @param B Weight specification used in the computation of the test statistic.
#'   One of \code{"identity"} (default), \code{"distribution"}, \code{"inverse diagonal"}, or \code{"inverse variance"} (not recommended).
#'   Alternatively, the user can supply a numeric matrix of dimension
#'   \code{length(mod$tau) x length(mod$tau)}. This argument is ignored when \code{full = TRUE}.
#' @param error.distr A character string specifying the assumed distribution of the
#'   regression errors, used only when \code{B = "distribution"}.
#'   Allowed values are \code{"normal"}, \code{"skew-normal"}, \code{"t"}.
#' @param error.par A named list of parameters associated with \code{error.distr}.
#'   Required elements depend on the chosen distribution:
#'   \itemize{
#'     \item For \code{"normal"}: \code{mean}, \code{sd}.
#'     \item For \code{"skew-normal"}: \code{xi}, \code{omega}, \code{alpha}.
#'     \item For \code{"t"}: \code{df}.
#'   }
#'
#' @return
#' A \code{data.frame} containing:
#' \itemize{
#'   \item \code{Quantiles.Set}: quantile levels
#'   \item \code{Statistic}: rank-score test statistic
#'   \item \code{p.value}: corresponding test-specific (unadjusted) \eqn{p}-value
#' }
#'
#' @author Angela Andreella
#'
#' @references
#' Koenker, R. (2005). \emph{Quantile Regression}. Cambridge University Press.
#'
#' Imhof, J. P. (1961). Computing the distribution of quadratic forms in normal
#' variables. \emph{Biometrika}, 48(3--4), 419--426.
#'
#' De Santis, R., Veselý, A., and Andreella, A. (2026). Inference on multiple
#' quantiles in regression models by a rank-score approach.
#' \emph{arXiv preprint} <doi:10.48550/arXiv.2511.07999>.
#'
#' @details
#' This procedure requires that the covariate of interest \code{X} is either numeric
#' or, if categorical, has at most two levels. Multilevel categorical covariates
#' are not supported and will trigger an error.
#'
#' If \code{full = TRUE}, \code{B} is ignored and the identity weight is used.
#' If \code{B} is user-defined (a matrix), it must be square, symmetric, positive
#' definite, and numerically invertible.
#'
#' @importFrom stats pgamma
#' @importFrom utils combn
#' @seealso \code{\link[quantreg]{rq}}, \code{\link[quasar]{waldTest}}
#' @export
#' @examples
#' set.seed(1234)
#' D <- simulateData(n = 100, gamma = 0.5, sigma.y = "1 + 2 * pmax(X, 0)")
#'
#' #Quantile regressions at different levels
#' tau <- c(0.1, 0.25, 0.5, 0.75, 0.9)
#' mod <- quantreg::rq(y ~ X + Z1, tau = tau, data=D)
#'
#' # Rank test
#' rankTest(mod, X = "X")

rankTest <- function(mod, X, tau = NULL, full = FALSE, h = NULL, alpha = 0.05, eps = c(1e-04,1e-04), B = "identity", error.distr = NULL, error.par = NULL){

  assert_intercept_present(mod = mod)
  assert_binary_categorical_X(mod = mod, X = X)

  if(is.null(tau)) tau <- mod$tau
  if (any(!tau %in% mod$tau)) stop("All values in tau must be among the quantiles used in mod.")


  res <- estimateCovariance(mod = mod, X = X, test = "rank", h = h, alpha = alpha)
  S <- res$S
  M <- res$M
  taus <- mod$tau


  # validate user-provided matrix B once (only if user supplies a matrix)
  if (is.matrix(B)) {

    if (nrow(B) != ncol(B)) {
      stop("Provided matrix 'B' must be square.", call. = FALSE)
    }

    if (nrow(B) != length(taus)) {
      stop(
        "Provided matrix 'B' must have dimension equal to the number of quantiles in mod$tau.",
        call. = FALSE
      )
    }

    # symmetry check
    if (!isTRUE(all.equal(B, t(B), tolerance = 1e-10))) {
      stop("Provided matrix 'B' must be symmetric.", call. = FALSE)
    }

    # full-rank / invertibility check (numerically stable)
    if (rcond(B) < .Machine$double.eps) {
      stop("Provided matrix 'B' must be full rank (numerically invertible).", call. = FALSE)
    }

    # positive definiteness check (no non-positive eigenvalues)
    ev <- eigen(B, symmetric = TRUE, only.values = TRUE)$values
    if (any(!is.finite(ev))) {
      stop("Provided matrix 'B' has non-finite eigenvalues.", call. = FALSE)
    }
    if (any(ev <= 0)) {
      stop("Provided matrix 'B' must be positive definite (all eigenvalues > 0).", call. = FALSE)
    }
  }


  tests <- unlist(lapply(1:length(taus),
                       combn,
                       x = taus,
                       simplify = FALSE),
                recursive = FALSE)

  pval  <- numeric(length(tests))
  tstat <- numeric(length(tests))
  set   <- character(length(tests))


  for (l in 1:length(tests)) {
    this_set <- unlist(tests[l])
    idx      <- which(taus %in% this_set)

    S_sub <- as.matrix(S[idx])
    M_sub <- as.matrix(M[idx, idx])


    if(B != "inverse variance"){
    if (identical(B, "identity")) {

      tstat[l] <- sum(S_sub**2)
      eigenvals <- eigen(M_sub, only.values = TRUE)$values

    } else if (identical(B, "inverse diagonal")) {

      tstat[l] <- sum((S_sub^2) * diag(1/M_sub))
      eigenvals <- eigen(M_sub %*% diag(diag(M_sub)**(-1), nrow = nrow(M_sub)), only.values = TRUE)$values

    } else if (identical(B, "distribution")) {

      if (is.null(error.distr)) {
        stop("Please specify the error distribution in 'error.distr'.", call. = FALSE)
      }

      B_vec <- build_B_from_dist(
        taus = this_set,
        error.distr = error.distr,
        error.par = error.par
      )

      tstat[l] <- sum((S_sub / B_vec)^2)
      eigenvals <- eigen(M_sub %*% (diag(B_vec^(-2), nrow = length(B_vec))), only.values = TRUE)$values

    } else if (is.matrix(B)) {

      B_sub <- as.matrix(B[idx, idx, drop = FALSE])

      tstat[l] <- t(S_sub) %*% solve(B_sub) %*% S_sub
      eigenvals <- eigen(M_sub %*% solve(B_sub), only.values = TRUE)$values

    } else {

      stop("'B' must be one of 'identity', 'distribution', 'inverse diagonal' or a matrix.", call. = FALSE)
    }


    if (any(Im(eigenvals) != 0)) {
      warning("Some eigenvalues are complex; using their real parts.", call. = FALSE)
    }
    eigenvals <- Re(eigenvals)
    eigenvals <- eigenvals[abs(eigenvals) > 1e-4]

    if (length(eigenvals) == 1) {
      pval[l] <- 1 - pgamma(tstat[l], shape = 1/2, scale = 2 * eigenvals)
    } else {
      pval[l] <- .pImhof(lams = eigenvals, x = tstat[l], eps = eps)
    }
  } else{
    tstat[l] <- (t(S_sub)%*%solve(M_sub)%*%S_sub)

    pval[l] <- 1-pchisq(tstat[l],df = length(S_sub))
  }
  }
  set = gsub(pattern = "c", replacement = "", x = paste0(tests))

  out <- data.frame(Quantiles.Set = set,
                    Statistic = tstat,
                    p.value = pval)

  if(full){
    pat <- paste0("(?<![0-9.])(", paste0(gsub("\\.", "\\\\.", tau), collapse="|"), ")(?![0-9])")
    out <- subset(out, grepl(pat, as.character(set), perl = TRUE))
  }else{
    pat <- paste0("^\\s*(", paste0(gsub("\\.", "\\\\.", tau), collapse="|"), ")\\s*$")
    out <- subset(out, grepl(pat, as.character(set)))
  }

  return(out)
}
