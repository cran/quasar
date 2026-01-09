#' @title Closed testing for quantile regression
#'
#' @description
#' Applies the closed testing procedure to strongly control the familywise error rate (FWER)
#' when testing the effect of a covariate of interest across multiple quantile regression models.
#'
#' @usage closedTesting(mod, X, tau = NULL, test = "rank-score", ...)
#'
#' @param mod An object of class \code{rqs} returned by
#'   \code{\link[quantreg]{rq}}, representing the fitted quantile regression models.
#' @param X A string indicating the covariate of interest.
#' @param tau A numeric vector of quantiles of interest used in \code{mod}.
#' If \code{NULL} (default), all quantiles from the \code{mod} object are considered.
#' @param test Character. Type of test to be used. Options are
#' \code{"rank-score"} and \code{"wald"}.
#' @param ... Additional arguments, see \code{\link[quasar]{rankTest}}, \code{\link[quasar]{waldTest}}.
#' @return
#' An object of class \code{quasar} containing:
#' \itemize{
#'   \item \code{Quantile}: quantile level
#'   \item \code{Coefficient}: estimated coefficient
#'   \item \code{Statistic}: test statistic
#'   \item \code{p.value}: raw \eqn{p}-value
#'   \item \code{p.value.adjusted}: adjusted \eqn{p}-value from the closed testing procedure
#' }
#'
#' @seealso
#' \code{\link[quantreg]{rq}}, \code{\link[quasar]{rankTest}}, \code{\link[quasar]{waldTest}}
#'
#' @references
#' Marcus, R., Eric, P., & Gabriel, K. R. (1976).
#' On closed testing procedures with special reference to ordered analysis of variance.
#' \emph{Biometrika}, 63(3), 655--660.
#'
#' Goeman, J. J., Hemerik, J., & Solari, A. (2021).
#' Only closed testing procedures are admissible for controlling false discovery proportions.
#' \emph{The Annals of Statistics}, 49(2), 1218--1238.
#'
#' @author Angela Andreella
#'
#' @details
#' This procedure requires that the covariate of interest \code{X} is either numeric
#' or, if categorical, has at most two levels. Multilevel categorical covariates
#' are not supported and will trigger an error.
#'
#' The weighting matrix used in the multivariate rank test is the identity matrix, i.e.,
#' it is currently the only one implemented explicitly as a shortcut within the closed testing procedure.
#'
#' @export
#'
#' @examples
#' # Simulate data
#' set.seed(1234)
#' D <- simulateData(n = 100, gamma = 0.5, sigma.y = "1 + 2 * pmax(X, 0)")
#'
#' # Quantile regressions at different levels
#' tau <- c(0.1, 0.25, 0.5, 0.75, 0.9)
#' mod <- quantreg::rq(y ~ X + Z1, tau = tau, data=D)
#'
#' # Closed testing
#' res <- closedTesting(mod, X = "X")
#' res
#'
#' # Summary and plot
#' summary(res, alpha = 0.1)
#' plot(res, alpha = 0.1, legend.position = "bottomright")

closedTesting <- function(mod, X, tau = NULL, test = "rank-score", ...){

  if (!inherits(mod, "rqs")) stop("mod must be an object of class rqs, typically returned by quantreg::rq().")

  if(is.null(tau)) tau <- mod$tau
  if (any(!tau %in% mod$tau)) stop("All values in tau must be among the quantiles used in mod.")

  test <- match.arg(test, c("rank-score", "wald"))

  parse_set <- function(s) {
    s <- gsub("[()\\s]", "", as.character(s))
    as.numeric(strsplit(s, ",")[[1]])
  }

  if (test == "rank-score") {
    res <- rankTest(mod = mod, X = X, tau = mod$tau, full = TRUE, ...)
  } else {
    res <- waldTest(mod = mod, X = X, tau = mod$tau, full = TRUE, ...)
  }

  out <- data.frame(p.value.adjusted = NA,
                    quantiles.set = NA)

  quantiles_single <- as.character(mod$tau)

  for(i in seq(length(quantiles_single))){
    idx <- vapply(res$Quantiles.Set, function(s) {
      vals <- parse_set(s)
      all(quantiles_single[i] %in% vals)
    }, logical(1))

    res_sub <- res[idx, ]
    out[i,] <- c(max(res_sub$p.value),
                 quantiles_single[i])
  }

  if(length(mod$tau)==1){
    out <- data.frame(Quantile = out$quantiles.set,
                      Coefficient = mod$coefficients[names(mod$model) == X],
                      Statistic = res$Statistic[1:length(mod$tau)],
                      p.value = round(res$p.value[1:length(mod$tau)],7),
                      p.value.adjusted = round(as.numeric(out$p.value.adjusted),7))
  }else{
    out <- data.frame(Quantile = out$quantiles.set,
                      Coefficient = mod$coefficients[names(mod$model) == X,],
                      Statistic = res$Statistic[1:length(mod$tau)],
                      p.value = round(res$p.value[1:length(mod$tau)],7),
                      p.value.adjusted = round(as.numeric(out$p.value.adjusted),7))

    rownames(out) <- NULL
    out <- subset(out, out$Quantile %in% tau)
  }

  class(out) <- c("quasar", class(out))
  return(out)
}
