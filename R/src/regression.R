ALPHA_XA <- function(data, response.vars, predictor.vars) {
  return(1)
}
ALPHA_XB <- function(data, response.vars, predictor.vars) {
  return(1)
}
ALPHA_AA <- function(data, response.vars, predictor.vars) {
  return(2)
}
ALPHA_AB <- function(data, response.vars, predictor.vars) {
  return(0.5)
}
B_A <- function(data, response.vars, predictor.vars) {
  return(1)
}
B_B <- function(data, response.vars, predictor.vars) {
  return(1)
}
COEF_MEAN_MEAN <- function(data, response.vars, predictor.vars) {
  preds <- obs.to.predictors(data[,predictor.vars,drop=F])
  return(array(0, c(ncol(preds)+1,length(response.vars))))
}
COEF_MEAN_COV <- function(data, response.vars, predictor.vars) {
  return(cov(data[,response.vars,drop=F]))
}
DATA_COV_SCALE <- function(data, response.vars, predictor.vars) {
  return(cov(data[,response.vars,drop=F]))
}
DATA_COV_DF <- function(data, response.vars, predictor.vars) {
  return(length(response.vars)+11)
}
COEF_COV_SCALE <- function(data, response.vars, predictor.vars) {
  preds <- obs.to.predictors(data[,predictor.vars,drop=F])
  phi <- diag(ncol(preds)+1)
  for (i in seq(2, length=ncol(preds))) {
    phi[i,i] <- 0.5
  }
  return(phi)
}
COEF_COV_DF <- function(data, response.vars, predictor.vars) {
  preds <- obs.to.predictors(data[,predictor.vars,drop=F])
  return(ncol(preds)+50)
}
COEF_COV_DF_LMB <- function(data, response.vars, predictor.vars) {
  preds <- obs.to.predictors(data[,predictor.vars,drop=F])
  return(ncol(preds)+0.5)
}

INIT_DPMLMB <- function(data, response.vars, predictor.vars) {
  preds <- obs.to.predictors(data[,predictor.vars,drop=F])
  d <- length(response.vars)
  h <- ncol(preds)+1
  return(list(
    A=array(0, c(h,d,1)),
    Sg=array(diag(d), c(d,d,1)),
    Omega=diag(h),
    M=array(0, c(h,d)),
    Z=as.integer(rep(0,nrow(data))),
    B=cbind(1, preds),
    x=0.5,
    al=1))
}

INIT_DPMLMLB <- function(data, response.vars, predictor.vars) {
  d <- length(response.vars)
  h <- 2
  return(list(
    A=array(0, c(h,d,1)),
    Sg=array(diag(d), c(d,d,1)),
    Omega=diag(h),
    M=array(0, c(h,d)),
    Z=as.integer(rep(0,nrow(data))),
    B=as.integer(rbinom(nrow(data), 1, 0.5)),
    x=0.5,
    al=1))
}

INIT_LMB <- function(data, response.vars, predictor.vars) {
  d <- length(response.vars)
  preds <- obs.to.predictors(data[,predictor.vars,drop=F])
  h <- ncol(preds)+1
  return(list(
    A=array(0, c(h,d)),
    Sg=diag(d),
    Omega=diag(h),
    B=cbind(1, preds)))
}

lmb <- function(response.vars, predictor.vars=NULL, class.var=NULL, data,
                      M=COEF_MEAN_MEAN(data, response.vars, predictor.vars),
                      Psi=DATA_COV_SCALE(data, response.vars, predictor.vars),
                      kappa=DATA_COV_DF(data, response.vars, predictor.vars),
                      Phi=COEF_COV_SCALE(data, response.vars, predictor.vars),
                      lambda=COEF_COV_DF_LMB(data, response.vars, predictor.vars),
                      init=INIT_LMB(data, response.vars, predictor.vars),
                      nsamp=500, nburnin=1200, lag=7, keep_chain=F) {
  x <- dataset.temp(response.vars, predictor.vars, class.var, data)
  preds <- obs.to.predictors(data[,predictor.vars,drop=F])
  hypers <- list(M=M,Psi=Psi,kappa=kappa,Phi=Phi,lambda=lambda)
  m <- jqgibbs(x, flgfm, "org/jqgibbs/models/MVRModel", NULL,
                       "org/jqgibbs/GenericSampler", hypers,
                       init, nsamp, nburnin, lag, NULL, F)
  m$obscol <- 1:ncol(preds)
  if (!keep_chain) {
    m$raw <- NULL
    gc()
  }
  return(m)
}

dpmlmb <- function(response.vars, predictor.vars=NULL, class.var=NULL, data,
                      W=COEF_MEAN_MEAN(data, response.vars, predictor.vars),
                      S=COEF_MEAN_COV(data, response.vars, predictor.vars),
                      Psi=DATA_COV_SCALE(data, response.vars, predictor.vars),
                      kappa=DATA_COV_DF(data, response.vars, predictor.vars),
                      Phi=COEF_COV_SCALE(data, response.vars, predictor.vars),
                      lambda=COEF_COV_DF(data, response.vars, predictor.vars),
                      xa=ALPHA_XA(data, response.vars, predictor.vars),
                      xb=ALPHA_XB(data, response.vars, predictor.vars),
                      ala=ALPHA_AA(data, response.vars, predictor.vars),
                      alb=ALPHA_AB(data, response.vars, predictor.vars),
                      init=INIT_DPMLMB(data, response.vars, predictor.vars),
                      nsamp=500, nburnin=1200, lag=7, keep_chain=F) {
  x <- dataset.temp(response.vars, predictor.vars, class.var, data)
  obscol <- match(predictor.vars, names(data))
  hypers <- list(W=W,S=S,Psi=Psi,kappa=kappa,Phi=Phi,lambda=lambda,
                 xa=xa,xb=xb,ala=ala,alb=alb,bshape1=1, bshape2=1)
  m <- jqgibbs(x, flgfa, "org/jqgibbs/models/FLGFDModel", NULL,
                       "org/jqgibbs/GenericSampler", hypers,
                       init, nsamp, nburnin, lag, NULL, F)
  m$obscol <- obscol
  if (!keep_chain) {
    m$raw <- NULL
    gc()
  }
  return(m)
}

dpmlmlb <- function(response.vars, predictor.vars=NULL, class.var=NULL, data,
                      W=COEF_MEAN_MEAN(data, response.vars, predictor.vars),
                      S=COEF_MEAN_COV(data, response.vars, predictor.vars),
                      Psi=DATA_COV_SCALE(data, response.vars, predictor.vars),
                      kappa=DATA_COV_DF(data, response.vars, predictor.vars),
                      Phi=COEF_COV_SCALE(data, response.vars, predictor.vars),
                      lambda=COEF_COV_DF(data, response.vars, predictor.vars),
                      xa=ALPHA_XA(data, response.vars, predictor.vars),
                      xb=ALPHA_XB(data, response.vars, predictor.vars),
                      ala=ALPHA_AA(data, response.vars, predictor.vars),
                      alb=ALPHA_AB(data, response.vars, predictor.vars),
                      bshape1=B_A(data, response.vars, predictor.vars),
                      bshape2=B_B(data, response.vars, predictor.vars),
                      init=INIT_DPMLMLB(data, response.vars, predictor.vars),
                      nsamp=500, nburnin=1200, lag=7, keep_chain=F) {
  x <- dataset.temp(response.vars, predictor.vars, class.var, data)
  obscol <- match(predictor.vars, names(data))
  hypers <- list(W=W,S=S,Psi=Psi,kappa=kappa,Phi=Phi,lambda=lambda,
                 xa=xa,xb=xb,ala=ala,alb=alb,be=bshape1,ga=bshape2)
  m <- jqgibbs(x, flgfe, "org/jqgibbs/models/FLGFEModel", NULL,
                       "org/jqgibbs/GenericSampler", hypers,
                       init, nsamp, nburnin, lag, NULL, F)
  m$obscol <- obscol
  if (!keep_chain) {
    m$raw <- NULL
    gc()
  }
  return(m)
}

dpmb <- function(response.vars, predictor.vars=NULL, class.var=NULL, data, ...) {
  return(dpmlmb(response.vars, predictor.vars, class.var, data, ...))
}

