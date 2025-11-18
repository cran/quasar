#' @importFrom Matrix bdiag
#' @importFrom quantreg rq
#' @importFrom quantreg ranks


estimateCovariance <- function(mod, X, test = "rank", h = NULL, alpha = 0.05){

  taus <- mod$tau
  density_list <- lapply(taus, function(tau) estimateDensity(mod = mod, tau, X = X, h = h, test = test, alpha = alpha))

  if(test == "rank"){

    f_quant_list <- lapply(density_list, function(x) 1/x)

    design <- mod$x
    formula <- make_h0_formula(mod = mod, X = X)
    design_h0 <- rq(formula, tau=taus, data = mod$model)$x

    #Proviamo con stima sotto H1 (24/10/2025)
  #  design_h0 <- mod$x

    x <- design[,(colnames(mod$model) %in% c(X))]

    hatz <-lapply(f_quant_list, function(w) design_h0%*%solve(t(design_h0)%*%diag(w)%*%design_h0)%*%t(design_h0)%*%diag(w)%*%x)

    a <- lapply(taus, function(w) ranks(rq(formula, tau=-1, data = mod$model), score="tau", tau=w))

    n <- nrow(design)

    s <- mapply(function(w, u) n^(-0.5)*t(x-w)%*%(u$ranks), hatz, a, SIMPLIFY = FALSE)

    Htaz <- do.call(cbind, hatz)

    S <- unlist(s)



    Z<-matrix(rep(x,length(taus)),n,length(taus))


    M<-t(Z-Htaz)%*%(Z-Htaz)/n

    ff<-Vectorize(function(x,y) {min(x,y)*(1-max(x,y))})
    M2<-outer(taus,taus,FUN=ff)
    M<-M*M2

  return(list(S = S, M = M))
  }else{

    design <- mod$x
    y <- mod$y
    z <- design[,!(colnames(design) %in% c("(Intercept)", X))]
    p <- ncol(design)
    n <- nrow(design)
    Tt <- length(taus)

    vcov <- matrix(0,Tt*p,Tt*p)

    J <- t(design)%*%(design)

    H <- sapply(seq(Tt), function (x) t(design)%*%(density_list[[x]]*design), simplify = FALSE)

    idx <- split(seq_len(p*Tt), rep(seq_len(Tt), each = p))

    Hinv <- lapply(H, solve)

    D <- bdiag(Hinv)
    K <- outer(taus, taus, function(a, b) pmin(a, b) - a*b)
    Vmid <- kronecker(K, J)
    vcov = as.matrix(D %*% Vmid %*% D)

    return(vcov = vcov)

  }
}
