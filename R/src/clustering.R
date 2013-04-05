###########################################################################
# Helper functions
############################################################################

l2norm <- function(v) {
  sqrt(sum(v^2))
}

glsap <- function(m) {
  Npreds <- ncol(m)
  if (ncol(m) < nrow(m)) {
    # This is the 'generalized' part...
    m <- cbind(m, matrix(max(m)+1, nrow(m), nrow(m)-ncol(m)))
  }
  o <- list()
  truth <- rownames(m)
  prediction <- colnames(m)
  s <- solve_LSAP(m)
  t <- as.numeric(s)
  t[t > Npreds] <- NA
  o$from.truth <- prediction[t]
  names(o$from.truth) <- truth
  o$from.prediction <- truth[order(t[!is.na(t)])]
  names(o$from.prediction) <- prediction[sort(t[!is.na(t)])]
  class(o) <- c("glsap")
  return(o)
}

precision <- function(pred, truth, pairwise=F) {
  if (pairwise) {
    pred <- samepairs(pred)
    truth <- samepairs(truth)
  }
  return(sum(pred & truth)/sum(pred))
}

recall <- function(pred, truth, pairwise=F) {
  if (pairwise) {
    pred <- samepairs(pred)
    truth <- samepairs(truth)
  }
  return(sum(pred & truth)/sum(truth))
}

fscore <- function(pred, truth, pairwise=F) {
  if (pairwise) {
    pred <- samepairs(pred)
    truth <- samepairs(truth)
  }
  p <- precision(pred, truth, pairwise=F)
  r <- recall(pred, truth, pairwise=F)
  f <- 2*p*r/(p+r)
  if (is.nan(f)) {
    return(0)
  } else {
    return(f)
  }
}

samepairs <- function(labels) {
  if (is.factor(labels)) {
    labels <- as.numeric(labels)
  }
  b <- outer(labels, labels, FUN="==")
  return(b[t(upper.tri(b))]) 
}

distmat.mlogdensity <- function(m, X, C) {
  d <- function(m, X) return(-density(m, X, log=T))
  return(distmat(m, X, C, distance=d, bias=T))
}

obs.to.predictor.col.indices <- function(obs, obscol=(1:ncol(obs))) {
  predictor.col.indices <- vector("list", length(obscol))
  names(predictor.col.indices) <- names(obs)[obscol]
  k <- 1
  for (o in names(obs)[obscol]) {
    if (is.factor(obs[[o]])) { 
      num_cols <- ncol(contrasts(obs[[o]]))
    } else {
      num_cols <- 1
    }
    predictor.col.indices[[o]] <- k:(k+num_cols-1)
    k <- k + num_cols
  }
  return(predictor.col.indices)
}

obs.to.predictors <- function(obs, obscol=(1:ncol(obs))) {
  predictor_cols <- obs.to.predictor.col.indices(obs, obscol)
  num_predictor_cols <- sum(as.numeric(sapply(predictor_cols, length)))
  predictors <- as.data.frame(matrix(0,nrow(obs),num_predictor_cols))
  for (o in names(obs)[obscol]) {
    v <- obs[[o]]
    if (is.factor(v)) {
      predictors[,predictor_cols[[o]]] <- contrasts(v)[v,]
    } else {
      predictors[,predictor_cols[[o]]] <- v
    }
  }
  return(predictors)
}

center.obs <- function(obs, obscol=(1:ncol(obs))) {
  predictors <- obs.to.predictors(obs, obscol)
  center <- rep(0, length(predictors))
  indices <- obs.to.predictor.col.indices(obs, obscol)
  for (o in names(obs)[obscol]) {
    v <- obs[[o]]
    if (is.factor(v)) {
      center[indices[[o]]] <- colMeans(contrasts(v))
    } else {
      center[indices[[o]]] <- mean(predictors[,indices[[o]]])
    }
  }
  return(center)
}

# Currently only handles factors
# Manipulate each variable separately, and let the remaining variables
# just be equal to their center value as given by center.obs
crucial.points <- function(obs, obscol=(1:ncol(obs))) {
  predictors <- obs.to.predictors(obs, obscol)
  indices <- obs.to.predictor.col.indices(obs, obscol)
  points <- t(predictors)
  points[,] <- center.obs(obs, obscol)
  num_points <- 0
  for (o in names(obs)[obscol]) {
    v <- obs[[o]]
    if (is.factor(v)) {
      num_points_o <- length(levels(v))
      points[indices[[o]],(num_points+1):(num_points+num_points_o)] <- t(contrasts(v)[levels(v),])
    } else {
      num_points_o <- 0
    }
    num_points <- num_points + num_points_o
  }
  return(t(points[,seq(1,length=num_points),drop=F]))
}

###########################################################################
# Generic functions
###########################################################################

evaluate <- function(m, X, C, ...) {
  UseMethod("evaluate")
}

#stats <- function(m, X, C, ...) {
#  UseMethod("stats")
#}
#
dataset <- function(d, ...) {
  UseMethod("dataset")
}

mixture <- function(...) {
  UseMethod("mixture")
}

component <- function(...) {
  UseMethod("component")
}

distmat <- function(m, X, C, ...) {
  UseMethod("distmat")
}

order_heur <- function(m) {
  UseMethod("order_heur")
}

########################################################################
# Dataset class
# Datasets are broken into three parts, by convention on the column
# labels: data (X), training-only data (T), and class labels (C)
# Training-only data is what we use to pass in the observable 'context'
# vectors, which don't get seen at test. 
########################################################################

dataset.temp <- function(response.vars, predictor.vars=NULL, class.var=NULL, data) {
  if (length(response.vars) == 0) {
    stop("No response variables selected")
  }
  ycols <- match(response.vars, names(data), nomatch=c())
  xcols <- match(predictor.vars, names(data), nomatch=c())
  ccol <- match(class.var, names(data), nomatch=c())
  if (length(ycols) != length(response.vars)) {
    stop("Unknown columns selected as response variables")
  }
  if (length(xcols) != length(predictor.vars)) {
    stop("Unknown columns selected as predictor variables")
  }
  if (length(class.var) > 0 & length(ccol) == 0) {
    stop("Unknown column selected as class variable")
  } else if (length(ccol) > 0) {
    ccol <- ccol[1]
  }
  data_new <- data[,c(ycols,xcols,ccol)]
  ycols_new <- match(response.vars, names(data_new), nomatch=c())
  xcols_new <- match(predictor.vars, names(data_new), nomatch=c())
  ccol_new <- match(class.var, names(data_new), nomatch=c())
  colnames(data_new)[ycols_new] <- paste("X", 1:length(ycols_new), sep='')
  if (length(xcols_new) > 0) {
    colnames(data_new)[xcols_new] <- paste("T", 1:length(xcols_new), sep='')
  }
  if (length(ccol_new) == 1) {
    colnames(data_new)[ccol_new] <- "C"
    data_new[["C"]] <- factor(data_new[["C"]])
  }
  return(dataset(data_new))
}

dataset.data.frame <- function(d) {
  # Create empty object
  ds <- list()
  # Skim identifying characters from column names
  ids <- substr(names(d), 1, 1)
  # Get data columns
  nda <- sort(names(d)[which(ids=="X")])
  ds$data <- d[,nda,drop=F]
  # Get training-only columns
  ntr <- sort(names(d)[which(ids=="T")])
  ds$tronly <- d[,ntr,drop=F]
  # Get class vector
  ncl <- sort(names(d)[which(ids=="C")])
  if (length(ncl) > 1) {
    ncl <- ncl[1]
  }
  ds$classes <- d[,ncl,drop=T]
  # Get secondary classification table
  ncl <- sort(names(d)[which(ids=="S")])
  if (length(ncl) != 0) {
    ds$secclasses <- d[,ncl,drop=F]
  }
  # Return object
  class(ds) <- "dataset"
  return(ds)
}

default.colnames <- function(ds) {
  colnames <- c()
  if (ncol(ds$data) > 0) {
    for (i in 1:ncol(ds$data)) {
      colnames <- c(colnames, paste("X",i,sep=''))
    }
  }
  if (ncol(ds$tronly) > 0) {
    for (i in 1:ncol(ds$tronly)) {
      colnames <- c(colnames, paste("T",i,sep=''))
    }
  }
  if (!identical(NULL, ds$classes)) {
    colnames <- c(colnames, "C")
  }
  if (!identical(NULL, ds$secclasses)) {
    for (i in 1:ncol(ds$secclasses)) {
      colnames <- c(colnames, paste("S",i,sep=''))
    }
  }
  return(colnames)
}

as.data.frame.dataset <- function(ds) {
  # Get default column names
  colnames <- default.colnames(ds)
  # Skim identifying characters from column names
  ids <- substr(colnames, 1, 1)
  # Construct empty data frame
  Ndata <- length(which(ids=="X"))
  Ntronly <- length(which(ids=="T"))
  Nsecs <- length(which(ids=="S"))
  d <- dataFrame(colClasses=
    c(rep("numeric", Ndata),
      rep("numeric", Ntronly),
      "numeric",
      rep("numeric", Nsecs)),
    nrow=nrow(ds$data)
  )
  names(d) <- colnames
  # Fill in data frame
  if (Ndata > 0) {
    d[,1:Ndata] <- ds$data
  }
  if (Ntronly > 0) {
    d[,(Ndata+1):(Ndata+Ntronly)] <- ds$tronly
  }
  if (!identical(NULL, ds$classes)) {
    d[,(Ndata+Ntronly+1)] <- ds$classes
  }
  if (!identical(NULL, ds$secclasses)) {
    d[,(Ndata+Ntronly+1+1):(Ndata+Ntronly+1+Nsecs)] <- ds$secclasses
  }
  # Return it
  return(d)
}

########################################################################
# Experiment class
########################################################################

expt <- function(name, setup) {
  e <- list(name=name, setup=setup)
  e$run <- function(...) { NULL }
  class(e) <- "expt"
  return(e)
}

########################################################################
# mclust Experiment subclass
########################################################################

mclustexpt <- function(name, setup) {
  e <- expt(name, setup)
  e$run <- function(data, obs) {
    settings <- e$setup(data, obs)
    j <- Mclust(data, settings$G, settings$modelNames, settings$prior)
    m <- mixture(j, data)
  }
  class(e) <- c(class(e), "mclustexpt")
  return(e)
}

###########################################################################
# Mixture model class
###########################################################################

mixture.default <- function(...) {
  l <- list()
  class(l) <- c("mixture")
  return(l)
}

#stats.mixture <- function(m) {
#  s <- rep(0, length(m$stats))
#  i <- 1
#  for (sn in names(m$stats)) {
#    s[i] <- m$stats[[sn]]()
#    i <- i + 1
#  }
#  names(s) <- names(m$stats)
#  return(s)
#}

evaluate.mixture <- function(m, X, C, metric=fscore, distmat=distmat.mlogdensity) {
  o <- list()
  # Get distance matrix
  o$distmat <- distmat(m, X, C)
  # Map classes in C to clusters in the model (and vice versa)
  g <- glsap(o$distmat)
  o$Ctom <- g$from.truth
  o$mtoC <- g$from.prediction
  # Remove the possibility of predicting a class not in Ctom (in case
  # there were too many clusters)
  n <- m
  n$active <- n$active %in% o$Ctom
  # Get predictions
  o$predictions <- o$mtoC[as.character(classify(m, X))]
  names(o$predictions) <- c()
  # Get performance
  o$performance <- list()
  for (cl in sort(unique(C))) {
    pred <- o$predictions == cl
    pred[is.na(pred)] <- FALSE
    truth <- C == cl
    o$performance[[cl]] <- metric(pred, truth)
  }
  pred <- o$predictions == C
  pred[is.na(pred)] <- FALSE
  truth <- C == C
  o$performance$overall <- metric(pred, truth)
  class(o) <- "evaluation"
  return(o)
}

mixture.Mclust <- function(o, data) {
  m <- mog(o, data)
  m$raw <- o
#  m$stats <- list(
#    BIC=function() {
#      m$raw$bic
#    },
#    K=function() {
#      length(m$active)
#    }
#  )
  return(m)
}

plot.mixture <- function(m,
                mixture_colors=default.mixture_colors,
                transform=default.transform,
                initplot=default.initplot, pch=20,
                comp_colors=NULL, ...) {
  # Initialize plot as directed
  if (!identical(initplot, FALSE)) {
    initplot(m, ...)
  }
  # Get color map
  colormap <- mixture_colors(m)
  # Plot points
  if (class(m)[1] %in% c("flgfd", "flgfa")) {
    z <- m$Z
  } else {
    z <- m$z
  }
  if (ncol(m$data$data) == 1) {
    for (k in unique(z)) {
      d <- m$data$data[z==k,1]
      if (length(d) > 1) {
        lines(density(m$data$data[z==k,1]), col=colormap[k])
      }
    }
  } else {
    pt_colors <- colormap[z]
    points(transform(m$data$data)[,1:2], col=pt_colors, pch=pch, ...)
  }
  # Plot components
  for (k in unique(z)) {
    if (identical(comp_colors,"mix")) {
      plot(component(m, k), m$data, col=colormap[k], transform=transform,
                            initplot=F, ...)
    } else if (!identical(comp_colors, NULL)) {
      plot(component(m, k), m$data, col=comp_colors, transform=transform,
                            initplot=F, ...)
    } else {
      plot(component(m, k), m$data, transform=transform,
                            initplot=F, ...)
    }
  }
}

order_heur.mixture <- function(m) {
  1:length(unique(m$Z))
}

###########################################################################
# Useful plotting code
###########################################################################
default.mixture_colors <- function(m) {
  cols <- rep(0, length(m$Z))
  used_cols <- unique(m$Z)
  cols[used_cols] <- order_heur(m)
  return(cols) 
}

default.transform <- function(d) {
  d[,1:2]
}

default.plotlim <- function(v, mar.sc=0.03) {
  vmin <- min(v)
  vmax <- max(v)
  vmar <- mar.sc*(abs(vmax-vmin))
  sgn.vmin <- vmin/abs(vmin)
  sgn.vmax <- vmax/abs(vmax)
  return(signif(c(vmin - vmar, vmax + vmar), 2))
}

default.initplot <- function(m, transform=default.transform,
                                plotlim=default.plotlim, ylim=c(0,2), ...) {
  plot.new()
  if (ncol(m$data$data) == 1) {
    xlim <- plotlim(m$data$data[,1])
  } else {
    xlim <- plotlim(transform(m$data$data)[,1])
    ylim <- plotlim(transform(m$data$data)[,2])
  }
  plot.window(xlim, ylim)
  axis(1)
  axis(2)
}
