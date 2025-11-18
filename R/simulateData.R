#' @title Simulate data
#' @description Simulates a main covariate \code{X},
#' a vector of additional covariates \code{Z}, and a response \code{y} drawn from
#' the chosen distribution.
#' @usage simulateData(n, beta = 0, gamma = 0, mu = 0, Sigma = NULL,
#'              sigma.y = 1, distribution = "normal", df = 5, seed = NULL)
#' @param n Integer. Number of observations.
#' @param beta Numeric scalar. Effect of \code{X}.
#' @param gamma Numeric vector. Effects of \code{Z} (length \code{p - 1}, where \code{p = ncol(Sigma)}).
#' @param mu Numeric scalar. Intercept.
#' @param Sigma Numeric \code{p x p} symmetric positive-definite covariance matrix for \code{(X, Z)}.
#'   The first column corresponds to \code{X}, the remaining columns to \code{Z1, Z2, ...}.
#' @param sigma.y Either a numeric scalar or a one-sided expression/string (e.g., \code{"0.3 * abs(X) + 0.1"})
#'   defining the scale of \code{y}.
#' @param distribution Character. One of \code{"normal"}, \code{"t"}, or \code{"exponential"}.
#'   This is the distribution of \code{y}.
#' @param df Numeric scalar > 0. Degrees of freedom for t-distribution.
#' @param seed Numeric scalar > 0. Seed for random number generator.
#' @author Angela Andreella
#' @return
#' A \code{data.frame} with columns \code{y}, \code{X}, and \code{Z1, ..., Zk}.
#'
#' @details
#' The response is generated as \code{y = mu + beta * X + Z \%*\% gamma + error}.
#' The error term can be drawn from a normal distribution, scaled Student-t with \code{df} degrees of freedom,
#' or a shifted exponential. Its standard deviation is defined by \code{sigma.y}:
#' if numeric, a fixed scale is used; if a character expression,
#' the scale can vary with \code{X} and/or \code{Z}.
#'
#' @examples
#' set.seed(1)
#' p <- 3
#' Sigma <- diag(p)
#'
#' # Normal
#' dat_n <- simulateData(n = 200, beta = 0.5, gamma = c(0.2,-0.1),
#'                       sigma.y = 0.5, distribution = "normal")
#'
#' # Student-t
#' dat_t0 <- simulateData(n = 200, beta = 0.5, gamma = c(0.2,-0.1),
#'                        sigma.y = 0.5, distribution = "t", df = 7)
#' # Exponential
#' dat_e <- simulateData(n = 200, beta = 0.5, gamma = c(0.2,-0.1),
#'                       sigma.y = "0.3 * abs(X) + 0.1", distribution = "exponential")
#'
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm rexp rt
#' @export

simulateData <- function(n, beta = 0, gamma = 0, mu = 0, Sigma = NULL, sigma.y = 1,
                         distribution = "normal", df = 5, seed = NULL) {

  if (!is.numeric(n) || length(n) != 1 || n <= 0) stop("n must be a positive integer.")
  n <- ceiling(n)

  if (!is.numeric(beta) || length(beta) != 1) stop("beta must be a numeric scalar.")
  if (!is.numeric(mu) || length(mu) != 1) stop("mu must be a numeric scalar.")
  if (!is.numeric(gamma)) stop("gamma must be a numeric vector.")

  p <- length(gamma) + 1

  if (is.null(Sigma)) Sigma <- diag(p)
  if (!is.matrix(Sigma) || ncol(Sigma) != nrow(Sigma)) stop("Sigma must be a square matrix.")
  if (!all(dim(Sigma) == p)) stop("Sigma dimensions must match length(gamma) + 1.")
  if (!isSymmetric(Sigma)) stop("`Sigma` must be symmetric.")

  distribution <- match.arg(tolower(distribution), c("normal", "t", "exponential"))

  if (is.numeric(sigma.y)) {
    if (any(sigma.y <= 0)) stop("All numeric values of sigma.y must be positive.")
  } else if (!is.character(sigma.y)) {
    stop("sigma.y must be either numeric or a character expression.")
  }

  if (!is.numeric(df) || length(df) != 1 || df <= 0) stop("df must be a positive number.")

  if (!is.null(seed)){
    if (!is.numeric(seed) || length(seed) != 1) stop("`seed` must be a numeric scalar.")
    set.seed(seed)
  }

  XZ <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(XZ) <- c("X", if (p > 1) paste0("Z", seq_len(p - 1)) else character(0))

  eta <- as.numeric(mu + beta * XZ[, 1] + if (p > 1) drop(XZ[, -1, drop = FALSE] %*% gamma) else 0)

  if(is.character(sigma.y)){
    env <- as.data.frame(XZ)
    sigma <- eval(str2lang(sigma.y), envir = env)
  }

  if (is.numeric(sigma.y)) {
    sigma <- rep(sigma.y, length.out=n)
  }

  if(distribution == "normal"){
    y <- rnorm(n = n, mean = eta, sd = sigma)
  }

  if (distribution == "t") {
    y <- eta + sigma * sqrt((df-2)/df) * rt(n = n, df = df)
  }

  if (distribution == "exponential") {
    eta_clip <- pmin(pmax(eta, -30), 30)
    rate <- 1 / (exp(eta_clip) * sigma)
    y <- rexp(n, rate = rate) - 1/rate
  }

  out <- as.data.frame(cbind(y, XZ))
  return(out)
}
