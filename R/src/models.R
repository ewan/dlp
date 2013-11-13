###########################################################################
# Useful functions
###########################################################################

all.subsets <- function(K, max) {
  v <- NULL
  for (k in 1:K) {
    if (is.null(v)) {
      mv <- 1
    } else {
      mv <- nrow(v)
    }
    v2 <- rbind(v,v)
    w <- rbind(matrix(0, mv, 1), matrix(1, mv, 1))
    v <- cbind(w, v2)
    s <- apply(v, 1, sum)
    v <- v[s<=max,,drop=F]
  }
  return(unique(v))
}

find.index.parms <- function(data, v, order, fixedparms, fixedlengths, dataparms,
                             fhparms, fhdims,
                             fdimparms, fdimdims, vdimparms, vdimdims,
                             datahparms, fdimhparms, vdimhparms, vhparms,
                             vhdims, h) {
  d <- ncol(data)
  N <- nrow(data)
  lengths <- fixedlengths
  names(lengths) <- fixedparms
  
  # Add data-length parms to lengths
  n <- names(lengths)
  for (parm in dataparms) {
    lengths <- c(lengths, nrow(data))
    n <- c(n, parm)
  }
  names(lengths) <- n

  # Add h-dependent data-length parms to lengths
  n <- names(lengths)
  for (parm in datahparms) {
    lengths <- c(lengths, h*N)
    n <- c(n, parm)
  }
  names(lengths) <- n

  # Add fixed-length h*h parms to lengths
  n <- names(lengths)
  for (parm in fhparms) {
    lengths <- c(lengths, h*h)
    n <- c(n, parm)
  }
  names(lengths) <- n

  # Add fixed-length h*dim parms to lengths
  n <- names(lengths)
  for (parm in fdimhparms) {
    lengths <- c(lengths, h*d)
    n <- c(n, parm)
  }
  names(lengths) <- n

  # Add fixed-length dim-dependent parms to lengths
  fdimlengths <- d^fdimdims
  n <- names(lengths)
  for (parm in fdimparms) {
    lengths <- c(lengths, fdimlengths[parm])
    n <- c(n, parm)
  }
  names(lengths) <- n

  # Work out K
  vdimtots <- c(d^vdimdims)
  vdimhtots <- rep(d*h, length(vdimhparms))
  vhtots <- c(h^vhdims)
  velems <- length(v) - sum(lengths)
  K <- velems/sum(sum(vdimtots),sum(vdimhtots),sum(vhtots))

  # Add variable-length h*dim parms to lengths
  n <- names(lengths)
  for (parm in vdimhparms) {
    lengths <- c(lengths, K*d*h)
    n <- c(n, parm)
  }
  names(lengths) <- n

  # Add variable-length dim-dependent parms to lengths
  vdimlengths <- K * vdimtots
  n <- names(lengths)
  for (parm in vdimparms) {
    lengths <- c(lengths, vdimlengths[parm])
    n <- c(n, parm)
  }
  names(lengths) <- n

  # Add variable-length h-dependent parms to lengths
  vhlengths <- K * vhtots
  n <- names(lengths)
  for (parm in vhparms) {
    lengths <- c(lengths, vhlengths[parm])
    n <- c(n, parm)
  }
  names(lengths) <- n
  # Get indexes
  index <- list()
  last.index <- 0
  for (parm in order) {
    end.index <- last.index + lengths[parm]
    index[[parm]] <- (last.index+1):end.index
    last.index <- end.index
  }

  return(index)
}

all.subsets <- function(K, max) {
  v <- NULL
  for (k in 1:K) {
    if (is.null(v)) {
      mv <- 1
    } else {
      mv <- nrow(v)
    }
    v2 <- rbind(v,v)
    w <- rbind(matrix(0, mv, 1), matrix(1, mv, 1))
    v <- cbind(w, v2)
    s <- apply(v, 1, sum)
    v <- v[s<=max,,drop=F]
  }
  return(unique(v))
}

###########################################################################
# Generic functions
###########################################################################

parm.index.finder <- function(m, data, v, order, ...) {
  UseMethod("parm.index.finder")
}

gaussian <- function(...) {
  UseMethod("gaussian")
}

bdev <- function(model, ...) {
  UseMethod("bdev")
}

pg <- function(...) {
  UseMethod("pg")
}

flgfm <- function(...) {
  UseMethod("flgfm")
}

lmod <- function(...) {
  UseMethod("lmod")
}

mog <- function(...) {
  UseMethod("mog")
}

###########################################################################
# Special numeric classes for parameters. These are used to define 'average.'
###########################################################################

numeric.mode <- function(v) {
  return(as.numeric(names(which.max(xtabs(~v)))))
}

numeric.incomplete <- function(x) {
  class(x) <- "numeric.incomplete"
  return(x)
}

average.numeric.incomplete <- function(x) {
  # FIXME - check dimensions somewhere so colMeans doesn't complain
  return(colMeans(x, na.rm=T))
}

discrete <- function(x) {
  class(x) <- "discrete"
  return(x)
}

average.discrete <- function(x) {
  mode <- c()
  for (i in 1:ncol(x)) {
    mode[i] <- numeric.mode(x[,i])
  }
  return(mode)
}

########################################################################
# Model classes and model function classes
########################################################################

flgfa.dataset <- function(ds) {
  o <- mixture()
  o$data <- ds$data
  o$A <- list()
  o$Sg <- list()
  o$B <- cbind(1, ds$tronly)
  clds <- sort(unique(as.numeric(ds$classes)))
  for (i in 1:length(clds)) {
    di <- o$data[as.numeric(ds$classes)==clds[i],,drop=F]
    bi <- o$B[as.numeric(ds$classes)==clds[i],,drop=F]
    b0i <- rep(T, nrow(di))
    for (j in 2:ncol(o$B)) {
      b0i <- b0i & !as.logical(bi[,j])
    }
    di0 <- di[b0i,]
    o$A[[i]] <- matrix(0, ncol(o$data), ncol(o$B))
    if (nrow(di0) > 0) {
      o$A[[i]][,1] <- colMeans(di0)
      o$Sg[[i]] <- cov(di0)
    } else {
      o$A[[i]][,1] <- colMeans(di)
      o$Sg[[i]] <- cov(di)
    }
    for (j in 2:ncol(o$B)) {
      dij <- di[bi[,j]==1,,drop=F]
      mn <- colMeans(dij)
      o$A[[i]][,j] <- mn - o$A[[i]][,1]
    }
  }
  o$Z <- as.numeric(ds$classes)
  o$active <- sort(unique(o$Z))
  class(o) <- c("flgfa", class(o))
  return(o)
}

flgfd.dataset <- function(ds) {
  o <- mixture()
  o$data <- ds
  o$A <- list()
  o$Sg <- list()
  o$B <- model.matrix(~ T1, data=ds$tronly) # FIXME
  o$obscol <- 1 # FIXME
  clds <- sort(unique(as.numeric(ds$classes)))
  for (i in 1:length(clds)) {
    di <- o$data$data[as.numeric(ds$classes)==clds[i],,drop=F]
    bi <- o$B[as.numeric(ds$classes)==clds[i],,drop=F]
    b0i <- rep(T, nrow(di))
    for (j in 2:ncol(o$B)) {
      b0i <- b0i & !as.logical(bi[,j])
    }
    di0 <- di[b0i,]
    o$A[[i]] <- matrix(0, ncol(o$data$data), ncol(o$B))
    if (nrow(di0) > 0) {
      o$A[[i]][,1] <- colMeans(di0)
      o$Sg[[i]] <- cov(di0)
    } else {
      o$A[[i]][,1] <- colMeans(di)
      o$Sg[[i]] <- cov(di)
    }
    for (j in 2:ncol(o$B)) {
      dij <- di[bi[,j]==1,,drop=F]
      mn <- colMeans(dij)
      o$A[[i]][,j] <- mn - o$A[[i]][,1]
    }
  }
  o$Z <- as.numeric(ds$classes)
  o$active <- sort(unique(o$Z))
  class(o) <- c("flgfa", class(o))
  return(o)
}

mog.dataset <- function(ds, predict_Z=F) {
  o <- mixture()
  o$data <- ds
  o$mu <- list()
  o$sg <- list()
  clds <- sort(unique(as.numeric(ds$classes)))
  for (c in clds) {
    di <- o$data$data[as.numeric(ds$classes)==c,,drop=F]
    o$mu[[c]] <- colMeans(di)
    o$sg[[c]] <- cov(di)
  }
  # Initially set actual assignments
  o$Z <- as.numeric(ds$classes)
  o$active <- sort(unique(o$Z))
  class(o) <- c("mog", class(o))
  # Replace assignments with model predictions
  if (predict_Z) {
    o$Z <- as.numeric(classify(o, o$data$data))
  }
  return(o)
}

mog.data.frame <- function(x, response_columns=NULL, class_column=NULL) {

  return(mog(dataset(x, response_columns, class_column)))
}

mixture.dataset <- function(ds, ...) {
  return(mog(ds, ...))
}

mog.Mclust <- function(m, data) {
  o <- mixture()
  if (!(missing(data))) {
    o$data <- data
    o$Z <- m$classification
    o$mu <- list()
    o$sg <- list()
    o$active <- sort(unique(o$Z))
    for (i in 1:length(o$active)) {
      o$mu[[i]] <- m$parameters$mean[,i]
      o$sg[[i]] <- m$parameters$variance$sigma[,,i]
    }
  }
  class(o) <- c("mog", class(o))
  return(o)
}

bdev.mog <- function(m, data=NULL) {
  if (is.null(data) | missing(data)) {
    data <- m$data
  }
  clusters <- as.character(m$active)
  freq <- summary(as.factor(m$Z))
  mds <- matrix(0, nrow(data), length(clusters))
  for (i in 1:length(clusters)) {
    mds[,i] <- density(gaussian(m, m$active[i]), data, log=T) + log(freq[clusters[i]])
  }
  loglik <- sum(apply(mds, 1, max))
  return(-2*loglik)
}

distmat.mog <- function(m, x, c, distance=mahal, bias=T) {
  classes <- as.character(sort(unique(c)))
  clusters <- as.character(m$active)
  freq <- summary(as.factor(m$Z))
  distmat <- matrix(0, length(classes), length(clusters))
  rownames(distmat) <- classes
  colnames(distmat) <- clusters
  for (cluster in clusters) {
    cluster.model <- gaussian(m, as.numeric(cluster))
    for (class in classes) {
      points <- x[c==class,]
      mds <- distance(cluster.model, points)
      if (bias) { # hmm... fixme
        mds <- mds + log(freq[cluster])
      }
      distmat[class,cluster] <- mean(mds)
    }
  }
  return(distmat)
}

classify.mog <- function(m, X) {
  clusters <- as.character(m$active)
  freq <- summary(as.factor(m$Z))
  mds <- matrix(0, nrow(X), length(clusters))
  for (i in 1:length(clusters)) {
    mds[,i] <- density(gaussian(m, m$active[i]), X, log=T) + log(freq[clusters[i]])
  }
  return(clusters[apply(mds, 1, which.max)])
}

component.mog <- function(m, i) {
  return(gaussian(m, i))
}

order_heur.mog <- function(m) {
  active <- unique(m$Z)
  l1 <- unlist(lapply(m$mu[active], sum))
  l2 <- unlist(lapply(m$mu[active], l2norm))
  return(order(l1+l2))
}

flgfe <- function(data, parms.to.elems, v, off.by.one=T) {
  o <- mixture()
  if (!(missing(data))) {
    d <- ncol(data$data)
    o$data <- data
    o$Z <- v[parms.to.elems[["Z"]]]
    if (off.by.one) {
      o$Z <- o$Z + 1
    }
    B <- v[parms.to.elems[["B"]]]
    o$B <- cbind(1, matrix(B, length(o$Z), 1, byrow=T))
    o$A <- list()
    o$Sg <- list()
    o$Omega <- list()
    A.raw <- v[parms.to.elems[["A"]]]
    Sg.raw <- v[parms.to.elems[["Sg"]]]
    h <- 2
    K <- length(A.raw[!is.na(A.raw)])/(d*h)
    for (i in 1:K) {
      offset.dh <- (i-1)*d*h
      offset.dd <- (i-1)*(d^2)
      offset.hh <- (i-1)*(h^2)
      o$A[[i]] <- matrix(A.raw[(offset.dh+1):(offset.dh+d*h)], d, h)
      o$Sg[[i]] <- matrix(Sg.raw[(offset.dd+1):(offset.dd+(d^2))], d, d)
    }
    o$Omega <- matrix(v[parms.to.elems[["Omega"]]], h, h)
    o$M <- matrix(v[parms.to.elems[["M"]]], d, h)
    o$x <- v[parms.to.elems[["x"]]]
    o$al <- v[parms.to.elems[["al"]]]
    o$active <- which(1:K %in% o$Z)
  }
  class(o) <- c("flgfe", "flgfd", class(o))
  return(o)
}


flgfa <- function(data, parms.to.elems, v, off.by.one=T) {
  o <- mixture()
  if (!(missing(data))) {
    d <- ncol(data$data)
    o$data <- data
    o$Z <- v[parms.to.elems[["Z"]]]
    if (off.by.one) {
      o$Z <- o$Z + 1
    }
    B <- v[parms.to.elems[["B"]]]
    h <- length(B)/length(o$Z)
    o$B <- matrix(B, length(o$Z), h, byrow=T)
    o$A <- list()
    o$Sg <- list()
    o$Omega <- list()
    A.raw <- v[parms.to.elems[["A"]]]
    Sg.raw <- v[parms.to.elems[["Sg"]]]
    Omega.raw <- v[parms.to.elems[["Omega"]]]
    K <- length(A.raw[!is.na(A.raw)])/(d*h)
    for (i in 1:K) {
      offset.dh <- (i-1)*d*h
      offset.dd <- (i-1)*(d^2)
      offset.hh <- (i-1)*(h^2)
      o$A[[i]] <- matrix(A.raw[(offset.dh+1):(offset.dh+d*h)], d, h)
      o$Sg[[i]] <- matrix(Sg.raw[(offset.dd+1):(offset.dd+(d^2))], d, d)
      o$Omega[[i]] <- matrix(Omega.raw[(offset.hh+1):(offset.hh+h^2)], h, h)
    }
    o$M <- matrix(v[parms.to.elems[["M"]]], d, h)
    o$x <- v[parms.to.elems[["x"]]]
    o$al <- v[parms.to.elems[["al"]]]
    o$active <- which(1:K %in% o$Z)
  }
  class(o) <- c("flgfa", class(o))
  return(o)
}

parm.index.finder.flgfa <- function(model, data, v, order, h) {
  dataparms <- c("Z")
  datahparms <- c("B")
  vdimparms <- c("Sg")
  vdimdims <- c(2)
  names(vdimdims) <- vdimparms
  vhparms <- c("Omega")
  vhdims <- c(2)
  names(vhdims) <- vhparms
  fdimparms <- c()
  fdimdims <- c()
  names(fdimdims) <- fdimparms
  fdimhparms <- c("M")
  vdimhparms <- c("A")
  lengths <- c(1, 1)
  names(lengths) <- c("al", "x")
  return(find.index.parms(data, v, order, names(lengths), lengths, dataparms, c(),
                          fdimparms, fdimdims, vdimparms, vdimdims, datahparms,
                          fdimhparms, vdimhparms, vhparms, vhdims, h))
}

bdev.flgfa <- function(m, data=NULL) {
  if (is.null(data) | missing(data)) {
    data <- m$data
  }
  clusters <- as.character(m$active)
  freq <- summary(as.factor(m$Z))
  mds <- matrix(0, nrow(data), length(clusters))
  for (i in 1:length(clusters)) {
    mds[,i] <- density(flgfm(m, m$active[i]), data, log=T) + log(freq[clusters[i]])
  }
  loglik <- sum(apply(mds, 1, max))
  return(-2*loglik)
}


wtfd2d <- function(o) {
	result <- matrix(o$rowVec()$value(), o$numRows(), o$numCols()) # byRow??
	return(result)
}

wtfd3d <- function(o) {
	a <- array(o$rowVec()$value(), c(o$numRows(), o$numCols(), o$size()))
	result <- lapply(apply(a, 3, list), function(m) m[[1]])
	return(result)
}

mlmp <- function(data, pteo, off.by.one=T) {
  o <- mixture()
  if (!(missing(data))) {
    d <- ncol(data$data)
    o$data <- data
    o$z <- as.numeric(pteo$get("z")$getNumericValue()$value())
    if (off.by.one) {
      o$z <- o$z + 1
    }
    o$A0 <- wtfd2d(pteo$get("A0")$getNumericValue())
    o$Omega <- wtfd2d(pteo$get("Omega")$getNumericValue())
    o$A <- wtfd3d(pteo$get("A")$getNumericValue())
    o$Sigma <- wtfd3d(pteo$get("Sigma")$getNumericValue())
    o$alpha <- pteo$get("alpha")$getNumericValue()$value()
    o$active <- unique(o$z)
  }
  class(o) <- c("mlmp", class(o))
  return(o)
}


#
#
#mlmp <- function(data, parms.to.elems, v, off.by.one=T) {
#  o <- mixture()
#  if (!(missing(data))) {
#    d <- ncol(data$data)
#    o$data <- data
#    o$z <- v[parms.to.elems[["z"]]]
#    if (off.by.one) {
#      o$z <- o$z + 1
#    }
#    h <- length(parms.to.elems[["A0"]])/d
#    o$A0 <- matrix(v[parms.to.elems[["A0"]]], d, h)
#    o$Omega <- matrix(v[parms.to.elems[["Omega"]]], h, h)
#    o$A <- list()
#    o$Sigma <- list()
#    A.raw <- v[parms.to.elems[["A"]]]
#    Sigma.raw <- v[parms.to.elems[["Sigma"]]]
#    K <- length(A.raw[!is.na(A.raw)])/(d*h)
#    for (i in 1:K) {
#      offset.dh <- (i-1)*d*h
#      offset.dd <- (i-1)*(d^2)
#      o$A[[i]] <- matrix(A.raw[(offset.dh+1):(offset.dh+d*h)], d, h)
#      o$Sigma[[i]] <- matrix(Sigma.raw[(offset.dd+1):(offset.dd+(d^2))], d, d)
#    }
#    o$alpha <- v[parms.to.elems[["alpha"]]]
#    o$active <- which(1:K %in% o$z)
#  }
#  class(o) <- c("mlmp", class(o))
#  return(o)
#}

parm.index.finder.mlmp <- function(model, data, v, order, h) {
  dataparms <- c("z")
  datahparms <- c()
  vdimparms <- c("Sigma")
  vdimdims <- c(2)
  names(vdimdims) <- vdimparms
  vhparms <- c()
  vhdims <- c()
  names(vhdims) <- vhparms
  fhparms <- c("Omega")
  fdimparms <- c()
  fdimdims <- c()
  names(fdimdims) <- fdimparms
  fdimhparms <- c("A0")
  vdimhparms <- c("A")
  lengths <- c(1)
  names(lengths) <- c("alpha")
  return(find.index.parms(data, v, order, names(lengths), lengths, dataparms,
                          fhparms,
                          fdimparms, fdimdims, vdimparms, vdimdims, datahparms,
                          fdimhparms, vdimhparms, vhparms, vhdims, h))
}

mlmvp <- function(data, parms.to.elems, v, off.by.one=T) {
  o <- mixture()
  if (!(missing(data))) {
    d <- ncol(data$data)
    o$data <- data
    o$z <- v[parms.to.elems[["z"]]]
    if (off.by.one) {
      o$z <- o$z + 1
    }
    h <- length(parms.to.elems[["A0"]])/d
    o$A0 <- matrix(v[parms.to.elems[["A0"]]], d, h)
    o$Omega <- matrix(v[parms.to.elems[["Omega"]]], h, h)
    o$A <- list()
    o$Sigma <- list()
    A.raw <- v[parms.to.elems[["A"]]]
    Sigma.raw <- v[parms.to.elems[["Sigma"]]]
    K <- length(A.raw[!is.na(A.raw)])/(d*h)
    for (i in 1:K) {
      offset.dh <- (i-1)*d*h
      offset.dd <- (i-1)*(d^2)
      o$A[[i]] <- matrix(A.raw[(offset.dh+1):(offset.dh+d*h)], d, h)
      o$Sigma[[i]] <- matrix(Sigma.raw[(offset.dd+1):(offset.dd+(d^2))], d, d)
    }
    o$alpha <- v[parms.to.elems[["alpha"]]]
    o$gamma <- v[parms.to.elems[["gamma"]]]
    o$active <- which(1:K %in% o$z)
  }
  class(o) <- c("mlmvp", class(o))
  return(o)
}

parm.index.finder.mlmvp <- function(model, data, v, order, h) {
  dataparms <- c("z")
  datahparms <- c()
  vdimparms <- c("Sigma")
  vdimdims <- c(2)
  names(vdimdims) <- vdimparms
  vhparms <- c()
  vhdims <- c()
  names(vhdims) <- vhparms
  fhparms <- c("Omega")
  fdimparms <- c()
  fdimdims <- c()
  names(fdimdims) <- fdimparms
  fdimhparms <- c("A0")
  vdimhparms <- c("A")
  lengths <- c(1, h)
  names(lengths) <- c("alpha", "gamma")
  return(find.index.parms(data, v, order, names(lengths), lengths, dataparms,
                          fhparms,
                          fdimparms, fdimdims, vdimparms, vdimdims, datahparms,
                          fdimhparms, vdimhparms, vhparms, vhdims, h))
}

order_heur.flgfa <- function(m) {
  active <- unique(m$Z)
  mu <- lapply(m$A[active], function(x) x[,1])
  l1 <- unlist(lapply(mu, sum))
  l2 <- unlist(lapply(mu, l2norm))
  return(order(l1+l2))
}

flgfm.default <- function(data, parms.to.elems, v, off.by.one=T) {
  o <- list()
  if (!(missing(data))) {
    d <- ncol(data$data)
    o$data <- data
    B <- v[parms.to.elems[["B"]]]
    h <- length(B)/nrow(data$data)
    o$B <- matrix(B, nrow(data$data), h, byrow=T)
    A <- v[parms.to.elems[["A"]]]
    o$A <- matrix(A, d, h)
    Sg <- v[parms.to.elems[["Sg"]]]
    o$Sg <- matrix(Sg, d, d)
    Omega <- v[parms.to.elems[["Omega"]]]
    o$Omega <- matrix(Omega, h, h)
  }
  class(o) <- c("flgfm")
  return(o)
}


flgfm.flgfa <- function(m, i) {
  o <- list()
  o$A <- m$A[[i]]
  o$Sg <- m$Sg[[i]]
  o$Omega <- m$Omega[[i]]
  o$B <- m$B
  o$d <- length(o$A) # FIXME wtf is this for
  o$h <- nrow(o$B)
  o$obscol <- m$obscol # FIXME also
  class(o) <- c("flgfm") 
  return(o)
}

component.flgfa <- function(m, i) {
  return(flgfm(m, i))
}

lmod.mlmp <- function(m, i) {
  o <- list()
  o$A <- m$A[[i]]
  o$h <- ncol(o$A)
  o$Sigma <- m$Sigma[[i]]
  o$d <- length(o$A)
  o$gamma <- rep(T, o$h)
  class(o) <- c("lmod") 
  return(o)
}

lmod.mlmvp <- function(m, i) {
  o <- list()
  o$gamma <- as.logical(m$gamma)
  o$A <- m$A[[i]][,o$gamma]
  o$h <- ncol(o$A)
  o$Sigma <- m$Sigma[[i]]
  o$d <- length(o$A)
  class(o) <- c("lmod") 
  return(o)
}

component.mlmp <- function(m, i) {
  return(lmod(m, i))
}

component.mlmvp <- function(m, i) {
  return(lmod(m, i))
}

density.flgfm <- function(m, X, log=T, maxc=4) {
  d <- c()
  if (ncol(m$B) > 1) {
    B <- all.subsets(ncol(m$B)-1, maxc)
    B <- cbind(rep(1, nrow(B)), B)
  } else {
    B <- 1
  }
  BA <- B%*%t(m$A)
  xMinusBA <- matrix(0, nrow(BA)*nrow(X), ncol(X))
  for (i in 1:ncol(X)) {
    xi <- rep(X[,i], nrow(BA))
    bai <- rep(0, nrow(BA)*nrow(X))
    for (j in 1:nrow(BA)) {
      bai[((j-1)*nrow(X)+1):(j*nrow(X))] <- BA[j,i]
    }
    xMinusBA[,i] <- xi - bai
  }
  d <- density(gaussian(rep(0, ncol(xMinusBA)), m$Sg), xMinusBA, log)
  dm <- matrix(d, nrow(X))
  return(apply(dm, 1, max))
}

plot.lmod <- function(m, data, transform=default.transform,
                          initplot=default.initplot,
                          lwd=1, lev=0.75, lty.cept=1,
                          lty.slope=2, npts1d=1000, old.school=F, ...) {
  if (!identical(initplot, FALSE)) {
    initplot(m)
  }
  if (old.school) {
#    h <- ncol(m$X) - 1
    h <- m$h-1
    center.x <- c(1, rep(0, h))
    center.pred <- predict(m, center.x)
    additional.points <- diag(h)
  } else {
    center.pred <- predict(m, c(1, center.obs(data$tronly))[m$gamma])
    additional.points <- crucial.points(data$tronly)
  }
  additional.preds <- predict(m, cbind(rep(1,nrow(additional.points)),
                                       additional.points)[,m$gamma])
#  if (nrow(m$Sigma) == 1) {
#    plotlim <- par("xaxp")[1:2]
#    plotpts <- seq(plotlim[1], plotlim[2], length=npts1d)
#    y <- dnorm(plotpts, center.pred, sqrt(m$Sigma))
#    lines(y ~ plotpts, lty=lty.cept, ...)
#    for (i in seq(1, length=nrow(additional.preds))) {
#      y <- dnorm(plotpts, additional.preds[i], sqrt(m$Sigma))
#      lines(y ~ plotpts, lty=lty.slope, ...)
#    }
#  } else {
    e <- ellipse(m$Sigma, centre = center.pred, lev = lev, ...)
    lines(transform(e[, 1:2]), lwd = lwd, lty = lty.cept, ...)
    for (i in seq(1, length=nrow(additional.preds))) {
      e <- ellipse(m$Sigma, centre = additional.preds[i,], lev = lev, lwd=lwd, ...)
      lines(transform(e[, 1:2]), lwd = lwd, lty = lty.slope, ...)
    }
#  }
}

predict.lmod <- function(m, X) {
  if (is.matrix(X)) {
    return(t(m$A %*% t(X)))
  } else {
    return(t(m$A %*% X))
  }
}

plot.flgfm <- function(m, data, axisstyle="normal", add=F,
                          xlim="auto", ylim="auto",
                          lwd=1, lev=0.75, lty.cept=1,
                          lty.slope=2, npts1d=1000, old.school=F, ...) {
  if (add==F) {
    plot.new()
    if (xlim == "auto") {
      xlim <- autolims(axisstyle, m)[[1]]
    }
    if (ylim == "auto") {
      ylim <- autolims(axisstyle, m)[[2]]
    }
    plot.window(xlim, ylim)
    axis(1)
    axis(2)
    title(xlab=xlab, ylab=ylab)
  }

  if (old.school) {
    h <- ncol(m$B) - 1
    center.x <- c(1, rep(0, h))
    center.pred <- predict(m, center.x)
    additional.points <- diag(h)
  } else {
    center.pred <- predict(m, c(1, center.obs(data$tronly)))
    additional.points <- crucial.points(data$tronly)
  }
  additional.preds <- predict(m, cbind(rep(1,nrow(additional.points)),
                                       additional.points))
  e <- ellipse(m$Sg, centre = center.pred, lev = lev, ...)
  lines(transform(e[,1:2], axisstyle), lwd = lwd, lty = lty.cept, ...)

  for (i in seq(1, length=nrow(additional.preds))) {
    e <- ellipse(m$Sg, centre = additional.preds[i,], lev = lev, lwd=lwd, ...)
    lines(transform(e[, 1:2], axisstyle), lwd = lwd, lty = lty.slope, ...)
  }
}

predict.flgfm <- function(m, B) {
  if (is.matrix(B)) {
    return(t(m$A %*% t(B)))
  } else {
    return(t(m$A %*% B))
  }
}

residuals.flgfm <- function(m, data) {
  pred <- predict(m, m$B)
  return(data$data-pred)
}

distmat.flgfa <- function(m, x, c, distance, bias=T) {
  classes <- as.character(sort(unique(c)))
  clusters <- as.character(m$active)
  freq <- summary(as.factor(m$Z))
  distmat <- matrix(0, length(classes), length(clusters))
  rownames(distmat) <- classes
  colnames(distmat) <- clusters
  for (cluster in clusters) {
    cluster.model <- flgfm(m, as.numeric(cluster))
    for (class in classes) {
      points <- x[c==class,]
      mds <- distance(cluster.model, points)
      if (bias) { # hmm... fixme
        mds <- mds + log(freq[cluster])
      }
      distmat[class,cluster] <- mean(mds)
    }
  }
  return(distmat)
}

classify.flgfm <- function(m, X, log=T, maxc=4) {
  d <- c()
  if (ncol(m$B) > 1) {
    B <- all.subsets(ncol(m$B)-1, maxc)
    B <- cbind(rep(1, nrow(B)), B)
  } else {
    B <- 1
  }
  B_ch <- apply(B[,2:ncol(B),drop=F], 1, function(x) paste(as.character(x), collapse='_'))
  BA <- B%*%t(m$A)
  xMinusBA <- matrix(0, nrow(BA)*nrow(X), ncol(X))
  for (i in 1:ncol(X)) {
    xi <- rep(X[,i], nrow(BA))
    bai <- rep(0, nrow(BA)*nrow(X))
    for (j in 1:nrow(BA)) {
      bai[((j-1)*nrow(X)+1):(j*nrow(X))] <- BA[j,i]
    }
    xMinusBA[,i] <- xi - bai
  }
  d <- density(gaussian(rep(0, ncol(xMinusBA)), m$Sg), xMinusBA, log)
  dm <- matrix(d, nrow(X))
  return(B_ch[apply(dm, 1, which.max)])
}


classify.flgfa <- function(m, X, predict_B=F, fixed_B=NULL, max=NULL) {
  clusters <- as.character(m$active)
  N <- nrow(X)
  h <- nrow(t(m$A[[m$active[1]]]))
  freq <- summary(as.factor(m$Z))
  if (!identical(max, NULL)) {
    if (max < length(clusters)) {
      freq <- freq[names(sort(freq, dec=T))[1:max]]
      clusters <- names(freq)
    }
  }
  mds <- matrix(0, nrow(X), length(clusters))
  if (!identical(fixed_B, NULL)) {
    for (i in 1:length(clusters)) {
      k <- as.numeric(clusters[i])
      corr <- as.matrix(fixed_B)[,2:h,drop=F] %*% t(m$A[[k]])[2:h,,drop=F]
      Xc <- X - corr
      mu <- t(m$A[[k]])[1,]
      sg <- m$Sg[[k]]
      mds[,i] <- density(gaussian(mu, sg), Xc, log=T) + log(freq[clusters[i]])
    }
    pred <- clusters[apply(mds, 1, which.max)]
    if (predict_B) {
      B_ch <- apply(as.matrix(fixed_B)[,2:h,drop=F], 1,
                    function(x) paste(as.character(x), collapse='_'))
      pred <- paste(pred, B_ch, sep='_')
    }
    return(pred)
  } else if (predict_B) {
    best_subs <- matrix('', nrow(X), length(clusters))
    for (i in 1:length(clusters)) {
      k <- as.numeric(clusters[i])
      Bi <- m$B[m$Z==k,,drop=F]
      B_ch <- apply(Bi[,2:ncol(Bi),drop=F], 1, function(x) paste(as.character(x), collapse='_'))
      sub_freq <- summary(as.factor(paste(i, B_ch, sep='_')))
      # Predict subcluster
      pred_B <- classify(flgfm(m, k), X) # N x h
      pred_sub <- paste(i, pred_B, sep='_')
      best_subs[,i] <- pred_sub
      # Correct data with respect to the appropriate subcluster
      pred_B_mat <- matrix(as.numeric(unlist(strsplit(pred_B, '_'))), N, h-1, byrow=T)
      corr <- pred_B_mat %*% t(m$A[[k]])[2:h,]
      Xc <- X - corr
      # Find the density of the data under the best subcluster assignment
      mu <- t(m$A[[k]])[1,]
      sg <- m$Sg[[k]]
      mds[,i] <- density(gaussian(mu, sg), Xc, log=T) + log(freq[clusters[i]])
                                                      + log(sub_freq[pred_sub])
    }
    cl_pred <- apply(mds, 1, which.max)
    preds <- rep("", N)
    for (i in 1:N) {
      preds[i] <- best_subs[i,cl_pred[i]]
    }
    return(preds)
  } else {
    for (i in 1:length(clusters)) {
      k <- as.numeric(clusters[i])
      mds[,i] <- density(flgfm(m, k), X, log=T) + log(freq[clusters[i]])
    }
    return(clusters[apply(mds, 1, which.max)])
  }
}


class(flgfa) <- "model.fn"
attr(flgfa, "types") <- list(A=numeric.incomplete, Sg=numeric.incomplete,
                        Omega=numeric.incomplete, Z=discrete, B=discrete,
                        x=numeric.incomplete, al=numeric.incomplete,
                        M=numeric.incomplete)
attr(flgfa, "params") <- c("A", "Sg", "Omega", "Z", "B", "x", "al", "M")

flgfd <- function(data, parms.to.elems, v, off.by.one=T) {
  o <- mixture()
  if (!(missing(data))) {
    d <- ncol(data$data)
    o$data <- data
    o$Z <- v[parms.to.elems[["Z"]]]
    if (off.by.one) {
      o$Z <- o$Z + 1
    }
    B <- v[parms.to.elems[["B"]]]
    h <- length(B)/length(o$Z)
    o$B <- matrix(B, length(o$Z), h, byrow=T)
    o$A <- list()
    o$Sg <- list()
    A.raw <- v[parms.to.elems[["A"]]]
    Sg.raw <- v[parms.to.elems[["Sg"]]]
    K <- length(A.raw[!is.na(A.raw)])/(d*h)
    for (i in 1:K) {
      offset.dh <- (i-1)*d*h
      offset.dd <- (i-1)*(d^2)
      o$A[[i]] <- matrix(A.raw[(offset.dh+1):(offset.dh+d*h)], d, h)
      o$Sg[[i]] <- matrix(Sg.raw[(offset.dd+1):(offset.dd+(d^2))], d, d)
    }
    o$Omega <- matrix(v[parms.to.elems[["Omega"]]], h, h)
    o$M <- matrix(v[parms.to.elems[["M"]]], d, h)
    o$x <- v[parms.to.elems[["x"]]]
    o$al <- v[parms.to.elems[["al"]]]
    o$active <- which(1:K %in% o$Z)
  }
  class(o) <- c("flgfd", class(o))
  return(o)
}

flgfm.flgfd <- function(m, i) {
  o <- list()
  o$A <- m$A[[i]]
  o$Sg <- m$Sg[[i]]
  o$B <- m$B
  o$d <- length(o$A)
  o$h <- nrow(o$B)
  class(o) <- c("flgfm") 
  return(o)
}

component.flgfd <- function(m, i) {
  return(flgfm(m, i))
}

classify.flgfd <- classify.flgfa

class(flgfd) <- "model.fn"
attr(flgfd, "types") <- list(A=numeric.incomplete, Sg=numeric.incomplete,
                        Z=discrete, B=discrete,
                        x=numeric.incomplete, al=numeric.incomplete,
                        Omega=numeric.incomplete, M=numeric.incomplete)
attr(flgfd, "params") <- c("A", "Sg", "Z", "B", "x", "al", "Omega", "M")


mpgha <- function(data, parms.to.elems, v, off.by.one=F) {
  o <- mixture()
  if (!missing(data)) {
    d <- ncol(data)
    o$data <- data
    o$mu <- list()
    o$be <- list()
    o$sg <- list()
    mu.raw <- v[parms.to.elems[["mu"]]]
    be.raw <- v[parms.to.elems[["be"]]]
    sg.raw <- v[parms.to.elems[["sg"]]]
    K <- length(mu.raw[!is.na(mu.raw)])/d
    for (i in 1:K) {
      offset.d <- (i-1)*d
      offset.dd <- (i-1)*(d^2)
      o$mu[[i]] <- mu.raw[(offset.d+1):(offset.d+d)]
      dim(o$mu[[i]]) <- c(1,d)
      o$be[[i]] <- be.raw[(offset.d+1):(offset.d+d)]
      dim(o$be[[i]]) <- c(1,d)
      o$sg[[i]] <- matrix(sg.raw[(offset.dd+1):(offset.dd+(d^2))], d, d)
    }
    o$Z <- v[parms.to.elems[["Z"]]]
    if (off.by.one) {
      o$Z <- o$Z + 1
    }
    o$B <- v[parms.to.elems[["B"]]]
    o$active <- which(1:length(o$mu) %in% o$Z)
    o$nmu <- v[parms.to.elems[["nmu"]]]
    o$nbe <- v[parms.to.elems[["nbe"]]]
    o$mu0 <- v[parms.to.elems[["mu0"]]]
    dim(o$mu0) <- c(1,ncol(o$data))
    o$be0 <- v[parms.to.elems[["be0"]]]
    dim(o$be0) <- c(1,ncol(o$data))
    o$x <- v[parms.to.elems[["x"]]]
    o$al <- v[parms.to.elems[["al"]]]
  }
  class(o) <- c("mpgha", class(o))
  return(o)
}
#
#bdev.mpgha <- function(m, data=NULL) {
#  if (is.null(data) | missing(data)) {
#    data <- m$data
#  }
#  clusters <- as.character(m$active)
#  freq <- summary(as.factor(m$Z))
#  mds <- matrix(0, nrow(data), length(clusters))
#  for (i in 1:length(clusters)) {
#    mds[,i] <- density(pg(m, m$active[i]), data, log=T) + log(freq[clusters[i]])
#  }
#  loglik <- sum(apply(mds, 1, max))
#  return(-2*loglik)
#}
#
#distmat.mpgha <- function(m, x, c, distance=mahal, bias=T) {
#  classes <- as.character(sort(unique(c)))
#  clusters <- as.character(m$active)
#  freq <- summary(as.factor(m$Z))
#  distmat <- matrix(0, length(classes), length(clusters))
#  rownames(distmat) <- classes
#  colnames(distmat) <- clusters
#  for (cluster in clusters) {
#    cluster.model <- pg(m, as.numeric(cluster))
#    for (class in classes) {
#      points <- x[c==class,]
#      mds <- distance(cluster.model, points)
#      if (bias) { # hmm... fixme
#        mds <- mds + log(freq[cluster])
#      }
#      distmat[class,cluster] <- mean(mds)
#    }
#  }
#  return(distmat)
#}
#
#classify.mpgha <- function(m, X) {
#  clusters <- as.character(m$active)
#  freq <- summary(as.factor(m$Z))
#  mds <- matrix(0, nrow(X), length(clusters))
#  for (i in 1:length(clusters)) {
#    mds[,i] <- density(pg(m, m$active[i]), X, log=T) + log(freq[clusters[i]])
#  }
#  return(clusters[apply(mds, 1, which.max)])
#}
#
#class(mpgha) <- "model.fn"
#attr(mpgha, "types") <- list(mu=numeric.incomplete, sg=numeric.incomplete,
#                        be=numeric.incomplete, Z=discrete, B=discrete,
#                        num=numeric.incomplete, nub=numeric.incomplete,
#                        mu0=numeric.incomplete, be0=numeric.incomplete,
#                        x=numeric.incomplete, al=numeric.incomplete)
#attr(mpgha, "params") <- c("mu", "sg", "be", "Z", "B", "num", "nub", "mu0", "be0", "x", "al")
#
#order_heur.mpgha <- order_heur.mog
#
#pg.mpgha <- function(m, i) {
#  o <- list()
#  o$mu <- m$mu[[i]]
#  o$be <- m$be[[i]]
#  o$sg <- m$sg[[i]]
#  o$pb <- sum(m$B)/length(m$B)
#  o$d <- length(o$mu)
#  class(o) <- c("pg") 
#  return(o)
#}
#
#density.pg <- function(m, X, log=T) {
#  if (is.null(dim(X))) {
#    if (length(X) == m$d+1) {
#      b <- as.numeric(X[length(X)])
#      return(density(gaussian(m$mu + b*m$be, m$sg), X[1:(length(X)-1)], log))
#    } else if (length(X) == m$d) {
#      m0 <- density(gaussian(m$mu, m$sg), X[1:(length(X))], log)
#      m1 <- density(gaussian(m$mu + m$be, m$sg), X[1:(length(X))], log)
#      return(max(m0,m1))
#    } else {
#      stop("Data has wrong number of dimensions")
#    }
#  } else {
#    return(apply(X, 1, function(x) {density.pg(m, x, log)}))
#  }
#}

#plot.pg <- function(m, initplot=default.initplot,
#                       ntps1d=1000,
#                       ...) {
#  if (!identical(initplot, FALSE)) {
#    initplot(m)
#  }
#  if (nrow(m$sg) == 1) {
#    plotlim <- par("xaxp")[1:2]
#    plotpts <- seq(plotlim[1], plotlim[2], length=npts1d)
#    y <- dnorm(plotpts, m$mu, sqrt(m$sg))
#    lines(y ~ plotpts, ...)
#    y <- dnorm(plotpts, m$mu+m$be, sqrt(m$sg))
#    lines(y ~ plotpts, ...)
#  } else {
#    e <- ellipse(m$sg, centre=m$mu, ...)
#    e <- ellipse(m$sg, centre=(m$mu+m$be), ...)
#    lines(transform(e[,1:2]), ...)
#    lines(transform(f[,1:2]), ...)
#  }
#}

gaussian.default <- function(mu, sg) {
  o <- list(mu=mu, sg=sg, d=length(mu))
  class(o) <- c("gaussian")
  return(o)
}

gaussian.mog <- function(m, i) {
  o <- list()
  o$mu <- m$mu[[i]]
  o$sg <- m$sg[[i]]
  o$d <- length(o$mu)
  class(o) <- c("gaussian") 
  return(o)  
}

density.gaussian <- function(m, X, log=T) {
  return(dmvnorm(X, m$mu, m$sg, log=log))
}

#plot.gaussian <- function (g, transform = default.transform, initplot = default.initplot, 
#    npts1d = 1000, ...) 
#{
#    if (!identical(initplot, FALSE)) {
#        initplot(g)
#    }
#    if (nrow(g$sg) == 1) {
#        plotlim <- par("xaxp")[1:2]
#        plotpts <- seq(plotlim[1], plotlim[2], length = npts1d)
#        y <- dnorm(plotpts, g$mu, sqrt(g$sg))
#        lines(y ~ plotpts, ...)
#    }
#    else {
#        em <- eigen(g$sg)
#        sm <- em$vectors
#        theta <- acos(g$sg[1,2]/(sqrt(g$sg[1,1]*g$sg[2,2])))
#        ellipse_shape(sqrt(qchisq(0.66, 2))*sqrt(g$sg[2,2]), 0,
#                      sqrt(qchisq(0.66, 2))*sqrt(g$sg[1,1]), 0,
#                      g$mu[c(2,1)], angle=theta, col=rgb(0.8,0.8,0.8,alpha=0.2))
#        e <- ellipse(g$sg, centre = g$mu, col='red', ...)
#        lines(transform(e[, 1:2]), lwd=0.05, ...)
#    }
#}


plot.gaussian <- function(g, add=F, npts1d=1000, xlim="auto", ylim="auto",
                             axisstyle="normal", ...) {
  if (add==F) { #FIXME
    plot.new()
    if (xlim == "auto") {
      xlim <- autolims(axisstyle, g)[[1]]
    }
    if (ylim == "auto") {
      ylim <- autolims(axisstyle, g)[[2]]
    }
    plot.window(xlim, ylim)
    axis(1)
    axis(2)
    title(xlab=xlab, ylab=ylab)
  }


  e <- ellipse(g$sg, centre=g$mu, ...) # FIXME
  lines(transform(e[,1:2], axisstyle), ...)
}

