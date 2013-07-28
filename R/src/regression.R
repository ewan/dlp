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
  if (is.null(predictor.vars)) {
    npreds <- 0
  } else {
    preds <- obs.to.predictors(data[,predictor.vars,drop=F])
    npreds <- ncol(preds)
  }
  return(array(0, c(npreds+1,length(response.vars))))
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
  if (is.null(predictor.vars)) {
    npreds <- 0
  } else {
    preds <- obs.to.predictors(data[,predictor.vars,drop=F])
    npreds <- ncol(preds)
  }
  phi <- diag(npreds+1)
  for (i in seq(2, length=npreds)) {
    phi[i,i] <- 0.5
  }
  return(phi)
}
COEF_COV_DF <- function(data, response.vars, predictor.vars) {
  if (is.null(predictor.vars)) {
    npreds <- 0
  } else {
    preds <- obs.to.predictors(data[,predictor.vars,drop=F])
    npreds <- ncol(preds)
  }
  return(npreds+50)
}
COEF_COV_DF_LMB <- function(data, response.vars, predictor.vars) {
  if (is.null(predictor.vars)) {
    npreds <- 0
  } else {
    preds <- obs.to.predictors(data[,predictor.vars,drop=F])
    npreds <- ncol(preds)
  }
  return(npreds+0.5)
}
TAU <- function(data, response.vars, predictor.vars) {
  return(0.5)
}

INIT_DPMLMVB <- function(data, response.vars, predictor.vars) {
  if (is.null(predictor.vars)) {
    npreds <- 0
    X <- data.frame(rep(1,nrow(data)))
  } else {
    preds <- obs.to.predictors(data[,predictor.vars,drop=F])
    npreds <- ncol(preds)
    X <- cbind(1, preds)
  }
  d <- length(response.vars)
  h <- npreds+1
  return(list(
    A=array(0, c(h,d,1)),
    Sigma=array(diag(d), c(d,d,1)),
    Omega=diag(h),
    A0=array(0, c(h,d)),
    z=as.integer(rep(0,nrow(data))),
    gamma=as.integer(rep(1,h)),
    X=X,
    alpha=1))
}

INIT_DPMLMB <- function(data, response.vars, predictor.vars) {
  if (is.null(predictor.vars)) {
    npreds <- 0
    X <- data.frame(rep(1,nrow(data)))
  } else {
    preds <- obs.to.predictors(data[,predictor.vars,drop=F])
    npreds <- ncol(preds)
    X <- cbind(1, preds)
  }
  d <- length(response.vars)
  h <- npreds+1
  return(list(
    A=array(0, c(h,d,1)),
    Sigma=array(diag(d), c(d,d,1)),
    Omega=diag(h),
    A0=array(0, c(h,d)),
    z=as.integer(rep(0,nrow(data))),
    X=X,
    alpha=1))
}

INIT_DPMLMLB <- function(data, response.vars, predictor.vars) {
  d <- length(response.vars)
  if (is.null(predictor.vars)) {
    h <- 1
  } else {
    h <- 2
  }
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

dpmlmvb <- function(response.vars, predictor.vars=NULL, class.var=NULL, data,
                      W=COEF_MEAN_MEAN(data, response.vars, predictor.vars),
                      S=COEF_MEAN_COV(data, response.vars, predictor.vars),
                      Psi=DATA_COV_SCALE(data, response.vars, predictor.vars),
                      kappa=DATA_COV_DF(data, response.vars, predictor.vars),
                      Phi=COEF_COV_SCALE(data, response.vars, predictor.vars),
                      lambda=COEF_COV_DF(data, response.vars, predictor.vars),
                      alpha_a=ALPHA_AA(data, response.vars, predictor.vars),
                      alpha_b=ALPHA_AB(data, response.vars, predictor.vars),
                      tau=TAU(data, response.vars, predictor.vars),
                      init=INIT_DPMLMVB(data, response.vars, predictor.vars),
                      nsamp=500, nburnin=1200, lag=7, keep_chain=F, 
                      T0z=1,Tfz=1e-3,T0g=1,Tfg=1, zlag=as.integer(1),
                      deadline=as.integer(nburnin),faster=F) {
  x <- dataset.temp(response.vars, predictor.vars, class.var, data)
  obscol <- match(predictor.vars, names(data))
  hypers <- list(W=W,S=S,Psi=Psi,kappa=kappa,Phi=Phi,lambda=lambda,
                 alpha_a=alpha_a,alpha_b=alpha_b,tau=tau,T0z=T0z,Tfz=Tfz,
                 T0g=T0g,Tfg=Tfg,deadline=deadline, zlag=zlag)
  if (faster == T) {
    model_class <- "org/jqgibbs/models/MLM_sample_params_varbsel_block_cache"
  } else {
    model_class <- "org/jqgibbs/models/MLM_sample_params_varbsel_block"
  }

  m <- jqgibbs(x, mlmvp, model_class, NULL, "org/jqgibbs/GenericSampler", hypers,
                       init, nsamp, nburnin, lag, NULL, F)
  m$obscol <- obscol
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
                      alpha_a=ALPHA_AA(data, response.vars, predictor.vars),
                      alpha_b=ALPHA_AB(data, response.vars, predictor.vars),
                      init=INIT_DPMLMB(data, response.vars, predictor.vars),
                      T0z=1,Tfz=1e-3,deadline=as.integer(nburnin),
                      nsamp=500, nburnin=1200, lag=7, keep_chain=F) {
  x <- dataset.temp(response.vars, predictor.vars, class.var, data)
  obscol <- match(predictor.vars, names(x$tronly))
  hypers <- list(W=W,S=S,Psi=Psi,kappa=kappa,Phi=Phi,lambda=lambda,T0z=T0z,Tfz=Tfz,
                 alpha_a=alpha_a,alpha_b=alpha_b, deadline=deadline)
  m <- jqgibbs(x, mlmp, "org/jqgibbs/models/MLM_sample_params", NULL,
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

