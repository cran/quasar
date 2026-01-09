
#' @importFrom stats as.formula
#' @importFrom stats formula
#' @importFrom stats terms
#' @importFrom pracma  quadgk

make_h0_formula <- function(mod, X) {
  f  <- formula(mod)
  tt <- terms(f)
  tl <- attr(tt, "term.labels")
  fac <- attr(tt, "factors")
  intc <- attr(tt, "intercept")

  if (length(tl) == 0L) {
    rhs <- if (intc == 1) "1" else "0"
    resp <- as.character(f)[2L]
    return(as.formula(paste(resp, "~", rhs)))
  }

  term_vars <- lapply(seq_along(tl), function(j) {
    rownames(fac)[fac[, j] != 0]
  })

  drop_idx <- which(vapply(term_vars, function(v) any(v %in% X), logical(1)))
  X_idx <- setdiff(seq_along(tl), drop_idx)

  rhs_terms <- tl[X_idx]
  rhs <- paste(c(if (intc == 1) "1" else "0", rhs_terms), collapse = " + ")
  rhs <- sub("^1 \\+ \\s*$", "1", rhs)
  rhs <- sub("^0 \\+ \\s*$", "0", rhs)

  resp <- as.character(f)[2L]
  as.formula(paste(resp, "~", if (nzchar(rhs)) rhs else if (intc == 1) "1" else "0"))
}

#' Assert that the fitted model includes an intercept
#'
#' @param mod An object of class "rq" or "rqs".
#' @keywords internal
assert_intercept_present <- function(mod) {
  has_int <- tryCatch({
    if(length(mod$tau)==1){
      rn <- names(stats::coef(mod))
    }else{
      rn <- rownames(stats::coef(mod))
    }
    "(Intercept)" %in% rn
  }, error = function(e) FALSE)

  if (!has_int) {
    stop("This method assumes a model with an intercept. ",
         "Please refit the quantile regression including an intercept (e.g., `y ~ 1 + ...`).",
         call. = FALSE)
  }
}


#' Assert that X is numeric or (if categorical) has at most two levels
#'
#' Works with rq/rqs even when model.frame or mod$x are not stored.
#' @param mod An "rq"/"rqs" object
#' @param X   Character scalar: name of the covariate of interest
#' @keywords internal

assert_binary_categorical_X <- function(mod, X) {
  stopifnot(is.character(X), length(X) == 1L)

  mf <- tryCatch(stats::model.frame(mod),
                 error = function(e) NULL)

  if (!is.null(mf) && X %in% names(mf)) {
    xv <- mf[[X]]
    check_logical <- ((is.factor(xv) && nlevels(xv) > 2L) | (is.character(xv) & length(unique(xv))>2L))
    levels_X <- ifelse(is.factor(xv), nlevels(xv), length(unique(xv)))
    if (check_logical) {
      stop(sprintf(
        "The covariate '%s' is a factor with %d levels. ",
        X, levels_X
      ),
      "This method only supports binary categorical covariates (<= 2 levels). ",
      "Please recode it before fitting the model.",
      call. = FALSE)
    }
    return(invisible(TRUE))
  }


}




#' Compute weights from error densities
#'
#' @param taus Numeric vector of quantile levels in (0,1).
#' @param error.distr Character. One of "normal", "skew-normal", "t".
#' @param error.par List of parameters for the chosen distribution.
#'
#' @keywords internal
#' @importFrom stats qnorm dnorm qt dt
#' @importFrom sn qsn dsn
build_B_from_dist <- function(taus, error.distr, error.par = list()) {
  error.distr <- match.arg(error.distr, c("normal", "skew-normal", "t"))

  # --- checks on taus ---
  if (!is.numeric(taus) || length(taus) == 0L) {
    stop("'taus' must be a non-empty numeric vector.", call. = FALSE)
  }
  if (anyNA(taus)) {
    stop("'taus' cannot contain NA values.", call. = FALSE)
  }
  if (any(taus <= 0 | taus >= 1)) {
    stop("'taus' must contain values strictly between 0 and 1.", call. = FALSE)
  }

  q_vals <- NULL
  f_vals <- NULL

  if (error.distr == "normal") {
    needed <- c("mean", "sd")
    missing <- needed[!needed %in% names(error.par)]
    if (length(missing) > 0L) {
      stop(
        sprintf("For 'normal' errors you must provide: %s.", paste(missing, collapse = ", ")),
        call. = FALSE
      )
    }

    mu <- error.par$mean
    sigma <- error.par$sd

    if (!is.numeric(mu) || length(mu) != 1L || is.na(mu)) {
      stop("'error.par$mean' must be a single numeric value.", call. = FALSE)
    }
    if (!is.numeric(sigma) || length(sigma) != 1L || is.na(sigma) || sigma <= 0) {
      stop("'error.par$sd' must be a single positive numeric value.", call. = FALSE)
    }

    q_vals <- stats::qnorm(taus, mean = mu, sd = sigma)
    f_vals <- stats::dnorm(q_vals, mean = mu, sd = sigma)

  } else if (error.distr == "skew-normal") {
    needed <- c("xi", "omega", "alpha")
    missing <- needed[!needed %in% names(error.par)]
    if (length(missing) > 0L) {
      stop(
        sprintf("For 'skew-normal' errors you must provide: %s.", paste(missing, collapse = ", ")),
        call. = FALSE
      )
    }

    xi    <- error.par$xi
    omega <- error.par$omega
    alpha <- error.par$alpha

    if (!is.numeric(omega) || length(omega) != 1L || is.na(omega) || omega <= 0) {
      stop("'error.par$omega' must be a single positive numeric value.", call. = FALSE)
    }

    q_vals <- sn::qsn(p = taus, xi = xi, omega = omega, alpha = alpha)
    f_vals <- sn::dsn(x = q_vals, xi = xi, omega = omega, alpha = alpha)

  } else { # error.distr == "t"
    if (is.null(error.par$df)) {
      stop("For 't' errors, you must provide 'df' in error.par.", call. = FALSE)
    }

    df <- error.par$df
    if (!is.numeric(df) || length(df) != 1L || is.na(df) || df <= 0) {
      stop("'error.par$df' must be a single positive numeric value.", call. = FALSE)
    }

    q_vals <- stats::qt(taus, df = df)
    f_vals <- stats::dt(q_vals, df = df)
  }

  if (any(!is.finite(f_vals))) {
    stop("Some densities are not finite. Check parameters and taus.", call. = FALSE)
  }
  if (any(f_vals <= 0)) {
    stop("All resulting densities must be positive.", call. = FALSE)
  }

  B <- (1 + f_vals)
  return(B)
}


#' The Imhof method for x=0 and central variables only.
#' @param lams eigenvalues
#' @param x observed stat test
#' @param eps tolerance vector (one for integrate, one for quadgk)
#' @importFrom stats integrate
#' @importFrom pracma quadgk
#' @importFrom methods is
.pImhof <- function(lams, eps = c(1e-04, 1e-04),  x) {
  lams <- lams[lams != 0]
  integrand <- function(u) {          # the Imhof integrand. Domain: 0...Inf
    theta <- 0.5 * colSums(atan(outer(lams,u)))  - 0.5 * x * u
    rho <- exp(colSums(0.25 * log(1 + outer(lams^2,u^2))))
    out <- ifelse(u==0, sum(lams)/2, sin(theta)/(u*rho))
    out
  }
  tr.integrand <- function(v) {       # Transformation of the integrand. Domain: 0...1
    K <- sum(abs(lams))/20            # Scaling constant of the transformation (to make it invariant to the scale of lams)
    0.5 + integrand(-log(1-v)/K) / (pi*K*(1-v))
  }
  rt1 <- max(sqrt(.Machine$double.eps), eps[1])
  rt2 <- max(sqrt(.Machine$double.eps), eps[2])
  res <- try(integrate(tr.integrand, 0, 1, rel.tol = rt1), silent=TRUE)
  if (is(res, "try-error")) {
    res <- try(quadgk(tr.integrand, 0, 1, tol = rt2), silent=TRUE)
  } else {
    out <- res$value
    return(out)
  }
  if (is(res, "try-error")) {
    out <- NA
  } else {
    out <- res
  }
  return(out)
}

