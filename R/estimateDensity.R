#' @importFrom quantreg rq
#' @importFrom quantreg bandwidth.rq
#' @importFrom quantreg rq.fit


estimateDensity <- function(mod, tau, X, test = "rank", h = NULL, alpha = 0.05){

  if(test == "rank"){
    #Proviamo con stima sotto H1 (24/10/2025) commento queste due righe
     formula <- make_h0_formula(mod = mod, X = X)
     mod<-rq(formula = formula, tau=tau, data = mod$model)
  }


  design <- mod$x
  y <- mod$y
  p <- ncol(design)
  n <- nrow(design)
  eps <- .Machine$double.eps^(1/2)
  resid <- mod$residuals

  pz <- sum(abs(resid) < eps)
  if(is.null(h) | identical(h, "Hall-Sheather")){
    h <- bandwidth.rq(p = tau, n = n, hs = TRUE, alpha = alpha)
  }
  if(identical(h, "Bofinger")){
    h <- bandwidth.rq(p = tau, n = n, hs = FALSE, alpha = alpha)
  }

  while ((tau - h < 0) || (tau + h > 1)) h <- h/2

  bhi <- rq.fit(design, y, tau = tau + h)$coef
  blo <- rq.fit(design, y, tau = tau - h)$coef
  dyhat <- design %*% (bhi - blo)

  f_quant_tau <- pmax(0.01, (2 * h)/(dyhat - eps))


    return(f_quant_tau)

}
