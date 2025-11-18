
#' @importFrom stats as.formula
#' @importFrom stats formula
#' @importFrom stats terms

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
