#' @title Rank-score test for quantile regression
#' @description
#' Performs the rank-score test for the covariate of interest
#' \code{X}, at the quantiles defined in \code{tau}, using a fitted quantile
#' regression model. The test evaluates the null hypothesis that the coefficient of \code{X} is equal to zero
#' against a two-sided alternative, at each specified quantile level.
#' Testing equality to a non-zero value is not yet implemented.
#'
#' @usage rankTest(mod, X, tau = NULL, full = FALSE, h = NULL, alpha = 0.05)
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
#' @param h A numeric value for the bandwidth.
#' @param alpha A numeric value used for bandwidth estimation.
#' Following Koenker (2005), it is typically set equal to the desired significance level.
#'
#' @return
#' A \code{data.frame} containing:
#' \itemize{
#'   \item \code{Quantiles.Set}: quantile levels
#'   \item \code{Statistic}: rank-score test statistic
#'   \item \code{p.value}: corresponding unadjusted \eqn{p}-value
#' }
#'
#' @author Angela Andreella
#'
#' @references
#' Koenker, R. (2005). \emph{Quantile Regression}. Cambridge University Press.
#'
#' @details
#' This procedure requires that the covariate of interest \code{X} is either numeric
#' or, if categorical, has at most two levels. Multilevel categorical covariates
#' are not supported and will trigger an error.
#'
#' @importFrom stats pchisq
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

rankTest <- function(mod, X, tau = NULL, full = FALSE, h = NULL, alpha = 0.05){

  assert_intercept_present(mod = mod)
  assert_binary_categorical_X(mod = mod, X = X)

  if(is.null(tau)) tau <- mod$tau
  if (any(!tau %in% mod$tau)) stop("All values in tau must be among the quantiles used in mod.")

  res <- estimateCovariance(mod = mod, X = X, test = "rank", h = h, alpha = alpha)
  S <- res$S
  M <- res$M
  taus <- mod$tau

  tests <- unlist(lapply(1:length(taus),
                       combn,
                       x = taus,
                       simplify = FALSE),
                recursive = FALSE)

  pval<-numeric(0)
  tstat <- numeric(0)

  for (l in 1:length(tests)) {
    S_sub <- as.matrix(S[which(taus %in% unlist(tests[l]))])

    M_sub <- as.matrix(M[which(taus %in% unlist(tests[l])),
                         which(taus %in% unlist(tests[l]))])

    tstat[l] <- (t(S_sub)%*%solve(M_sub)%*%S_sub)

    pval[l] <- 1-pchisq(tstat[l],df = length(S_sub))
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
