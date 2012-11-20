###############################################################################
# qgibbs.R
# Useful functions for Gibbs sampling
###############################################################################
# 14/12/2010: Collected the old code to read and write chains in my format
#             plus the code to control my Java sampler here in this file
###############################################################################
# Copyright Ewan Dunbar
# University of Maryland
##############################################################################


##############################################################################
# Constants
# FIXME: these should be more customizable
##############################################################################

QG.MODEL.FNS <- list(
  flgff=flgfa,
  flgfc=flgfa,
  flgfa=flgfa,
  flgfe=flgfe,
  mpg=mpgha
)

JQG.MODEL.CLASSES <- list(
  flgff="org/jqgibbs/models/FLGFFModel",
  flgfc="org/jqgibbs/models/FLGFCModel",
  flgfa="org/jqgibbs/models/FLGFAModel",
  mpg="org/jqgibbs/models/MpgFModel"
)

##############################################################################
# Generic functions
##############################################################################

bind <- function(x, y) {
  UseMethod("bind")
}

DIC <- function(x) {
  UseMethod("DIC")
}

acor <- function(x, ...) {
  UseMethod("acor")
}

##############################################################################
# Run jqgibbs models
##############################################################################

jqgibbs <- function(data, model.fn, model.class, name, samplerClass, hypers, init,
                    nsamples, nburnin, lag, savetime=NULL, finalsave=F) {
  if (is.null(names(hypers)) | is.null(names(init))) {
    stop("Both hypers and init must have names")
  }
  # Initialize JVM
  .jinit()
  .jaddClassPath(paste(.path.package("dlp"), "java", sep="/"))
  # Prepare the Java parameters
  dm <- as.matrix(data$data)
  storage.mode(dm) <- "double"
  dataobj <- jqgobj(dm)
  hyperobj <- jqgobj(hypers)
  initobj <- jqgobj(init)
  # Prepare the sampler
  modelobj <- .jnew(model.class, hyperobj, initobj, dataobj$numCols())
  modelobj <- .jcast(modelobj, "org/jqgibbs/Model")
  samplerobj <- .jnew(samplerClass, modelobj, dataobj)
  # Run the sampler, saving about every savetime
  lastSaveTime <- proc.time()["elapsed"] 
  chobj <- .jnew("org/jqgibbs/Chain")
  clean <- TRUE
  i <- 0
  realiters <- 0
  totaliters <- nburnin + lag*nsamples
  offiters <- 0
  burnedIn <- !(nburnin > 0)
  lagIter <- FALSE
  pbar <- txtProgressBar(min = 0, max = totaliters, style = 3)
  while (i < nsamples) {
    chlinkobj <- samplerobj$variateFast()
    currTime <- proc.time()["elapsed"]
    if (!(is.null(savetime))) {
      if ((i >= 1) & (currTime-lastSaveTime > savetime)) {
        # Save to disk
        ch <- qchain(chobj, dm, model.fn, model.class)
        save.qchain(ch, name)
        # Remove chain object
        rm(ch)
        gc()
        # Reset save time
        lastSaveTime <- currTime
        clean <- TRUE
      }
    }
    # Increment burnin, lag, and sample counters
    if (burnedIn) {
      if (lag > 0 & lagIter) {
        offiters <- offiters + 1
        if (offiters > lag) {
          stop("BUG: should never count more than lag steps between samples")
        }
        if (offiters == lag) {
          lagIter <- FALSE
        }
      } else {
        chobj$addLink(chlinkobj)
        clean <- FALSE
        i <- i + 1
        offiters <- 0
        lagIter <- TRUE
      }
    } else {
      offiters <- offiters + 1
      if (offiters > nburnin) {
        stop("BUG: should never count more than burnin before sample")
      }
      if (offiters == nburnin) {
        burnedIn <- TRUE
        offiters <- 0
      }
    }
    # Update progress bar
    realiters <- realiters + 1
    setTxtProgressBar(pbar, realiters)
  }
  # Finish off that progress bar
  write("\n", file="")
  # Save the chain a final time
  if (finalsave) {
    ch <- qchain(chobj, dm, model.fn, model.class)
    save.qchain(ch, name)
    rm(ch)
  }
  # Return the mixture object
  m <- mixture(chobj, data, model.fn, model.class)
  rm(chobj)
  gc()
  return(m)
}

############################################################################
# Generic functions
############################################################################

jqgobj <- function(x, ...) {
  UseMethod("jqgobj")
}

jqgclass <- function(x, ...) {
  UseMethod("jqgclass")
}

qchain <- function(x, data, model.fn, ...) {
  UseMethod("qchain")
}

average <- function(x) {
  UseMethod("average")
}

##############################################################################
# Create jqgibbs objects
##############################################################################

jqgobj.integer <- function(x) {
  if (length(x) == 1) {
    return(.jnew("org/jqgibbs/mathstat/Integer0D", x))
  } else {
    return(.jnew("org/jqgibbs/mathstat/Integer1D", x))
  }
}

jqgobj.numeric <- function(x) {
  if (is.null(dim(x)) | length(dim(x)) == 1) {
    if (length(x) == 1) {
      return(.jnew("org/jqgibbs/mathstat/Double0D", x))
    } else {
      return(.jnew("org/jqgibbs/mathstat/Double1D", x))
    }
  }
}

jqgobj.data.frame <- function(x) {
  .jinit()
  return(jqgobj(as.matrix(x)))
}

jqgobj.array <- function(x) {
  .jinit()
  if (length(dim(x)) == 3) {
    y <- array(0, dim(x)[c(2,1,3)])
    for (i in 1:dim(x)[3]) {
      y[,,i] <- t(x[,,i])
    }
    l <- as.integer(dim(x)[3])
    m <- as.integer(dim(x)[1])
    n <- as.integer(dim(x)[2])
    v <- as.vector(y)
    if (length(v) > 1) {
      if (is.integer(v)) {
        obj <- .jnew("org/jqgibbs/mathstat/Integer3D", l, m, n, v)
      } else {
        obj <- .jnew("org/jqgibbs/mathstat/Double3D", l, m, n, v)
      }
    } else {
      obj <- jqgobj(v) 
      obj <- obj$sequence()$sequence()$sequence()
    }
    return(obj)
  } else {
    # FIXME
    # It's complicated. Use extract with an index array.
  }
}

jqgobj.matrix <- function(x) {
  .jinit()
  m <- as.integer(nrow(x))
  n <- as.integer(ncol(x))
  v <- as.vector(t(x))
  if (length(v) > 1) {
    if (is.integer(v)) {
      obj <- .jnew("org/jqgibbs/mathstat/Integer2D", m, n, v)
    } else {
      obj <- .jnew("org/jqgibbs/mathstat/Double2D", m, n, v)
    }
  } else {
    obj <- jqgobj(v) 
    obj <- obj$sequence()$sequence()
  }
  return(obj)
}

jqgobj.list <- function(x) {
  .jinit()
  if (is.null(names(x))) {
    obj <- .jnew("java/util/LinkedList")
    for (y in x) {
      obj$add(jqgobj(y))
    }
  } else {
    obj <- .jnew("java/util/HashMap")
    for (i in 1:length(x)) {
      n <- names(x)[i]
      y <- x[[i]]
      obj$put(n, jqgobj(y))
    }
    obj <- .jcast(obj, "java/util/Map")
  }
  return(obj)
}

jqgobj.model.fn <- function(x, hyperobj, initobj, dataobj) {
  .jinit()
  model.class <- jqgclass(x)
  obj <- .jnew(model.class, hyperobj, initobj, dataobj$numCols())
  obj <- .jcast(obj, "org/jqgibbs/Model")
  return(obj)
}

##############################################################################
# Get R objects from jqgobj objects
##############################################################################

# FIXME
# There's actually no class called jqgobj! So I'm calling these directly
# until I figure out how to extend the S4 jobjRef objects seamlessly.

as.list.jqgobj <- function(x, ...) {
  if (.jinstanceof(x, "java/util/Map")) {
    entries <- as.list(x$entrySet())
    l <- vector("list", length(entries))
    n <- vector("character", length(entries))
    for (i in 1:length(entries)) {
      entry <- entries[[i]]
      key <- entry$getKey()
      if (class(key) == "jobjRef") {
        key <- key$toString()
      }
      value <- entry$getValue()
      if (class(value) == "jobjRef") {
        value <- as.numeric.jqgobj(value)
      }
      n[i] <- key
      l[[i]] <- value
    }
    names(l) <- n
    return(l)
  } else {
    # Here you would remove the outer layer of structure - FIXME
    return(as.list(x))
  }
}

as.numeric.jqgobj <- function(x, ...) {
  if (.jinstanceof(x, "org/jqgibbs/mathstat/Numeric")) {
    v <- x$value()
    if (class(v) == "jobjRef") {
      return(as.numeric.jqgobj(v))
    }
    return(as.numeric(v))
  } else {
    # Here you would remove the outer layer of structure - FIXME
    return(as.numeric(x))
  }
}

##############################################################################
# Retrieve jqgibbs class names
##############################################################################

jqgclass.model.fn <- function(x) {
  if (!(list(x) %in% QG.MODEL.FNS)) {
    stop("Unknown model function")
  }
  i <- match(list(x), QG.MODEL.FNS)
  s <- names(QG.MODEL.FNS)[i]
  if (!(s %in% names(JQG.MODEL.CLASSES))) {
    stop("Unsupported model function")
  }
  return(JQG.MODEL.CLASSES[s])
}

##############################################################################
# Chain class
##############################################################################

qchain.old.style <- function(x, data, model.fn, ...) {
  ch <- list()
  ch$data <- data
  ch$model.fn <- model.fn
  # Remove the old-style params which aren't model parms, but keep order
  params <- attr(x, "params")
  relevant <- params %in% attr(model.fn, "params")
  relevant.params <- params[relevant]
  # Check the converse 
  not.missing.parms <- as.logical(prod(attr(model.fn, "params") %in% params))
  if (!not.missing.parms) {
    stop("Not all model parameters are in the old-style chain")
  }
  # Work out which columns should correspond to which parameters
  n.rows <- length(x)
  ch$columns <- list()
  last.column <- 0
  for (i in 1:length(params)) {
    current.last <- last.column + length(unlist(x[[n.rows]][[i]]))
    p <- params[i]
    if (p %in% relevant.params) {
      ch$columns[[p]] <- (last.column+1):current.last
      last.column <- current.last
    }
  }
  n.cols <- last.column
  # Fill in matrix
  ch$m <- matrix(nrow=n.rows, ncol=last.column)
  for (i in 1:length(x)) {
    for (p in relevant.params) {
      j <- match(p, params)
      n.elements <- length(unlist(x[[i]][[j]]))
      columns.par <- (ch$columns[[p]][1]):(ch$columns[[p]][1]+n.elements-1)
      ch$m[i,columns.par] <- unlist(x[[i]][[j]])
    }
  }
  # Return the chain
  ch$off.by.one <- F
  ch$dim <- ncol(ch$data)
  class(ch) <- "qchain"
  return(ch)
}

qchain.character <- function(x, data, model.fn, parm.order, ...) { # For dumpfiles
  ch <- NULL
  lines <- readLines(x)
  links <- lapply(strsplit(lines, " "), as.numeric)
  for (link in links) {
    ch.curr <- list()
    class(ch.curr) <- "qchain"
    ch.curr$data <- data
    ch.curr$model.fn <- model.fn
    ch.curr$columns <- parm.index.finder(model.fn(), data, link, parm.order, ...)
    ch.curr$m <- matrix(link, 1)
    if (is.null(ch)) {
      ch <- ch.curr
    } else {
      ch <- bind(ch, ch.curr)
    }
  }
  # Return the chain
  ch$off.by.one <- T # FIXME - not good enough for the long chain summary
  ch$dim <- ncol(ch$data)
  return(ch)
}

# FIXME - change this when you figure out how to create the jqgobj
# class
qchain.jobjRef <- function(x, data, model.fn, jqgmodel, ...) {
  ch <- NULL
  ch$data <- data
  ch$model.fn <- model.fn
  ch$jqgmodel <- jqgmodel
  ch$m <- x$getFlatDouble2D()$value()
  ch$columns <- as.list.jqgobj(x$getOrderMap())
  for (n in names(ch$columns)) {
    ch$columns[[n]] <- ch$columns[[n]] + 1
  }
  ch$off.by.one <- T # FIXME - not good enough for the long chain summary
  ch$dim <- ncol(ch$data)
  class(ch) <- "qchain"
  return(ch)
}

"[[.qchain" <- function(ch, i) {
   if (length(i) != 1) {
    stop("Cannot select multiple model objects from chain")
   }
   return(ch$model.fn(ch$data, ch$columns, ch$m[i,], ch$off.by.one))
}

"[.qchain" <- function(ch, i) {
  chp <- ch
  chp$m <- chp$m[i,,drop=F]
  return(chp)
}

length.qchain <- function(c) {
  return(nrow(c$m))
}

bind.qchain <- function(c1, c2) {
  if (class(c1$model.fn) != class(c2$model.fn)) {
    stop("Attempt to bind incompatible chains")
  }
  equal <- F
  if (ncol(c1$m) > ncol(c2$m)) {
    bigger <- c1
    smaller <- c2
  } else if (ncol(c1$m) < ncol(c2$m)) {
    bigger <- c2
    smaller <- c1
  } else {
    equal <- T
  }
  if (!equal) {
    new.smaller.m <- matrix(nrow=length(smaller), ncol=ncol(bigger$m))
    columns <- bigger$columns
    params <- names(columns)
    for (p in params) {
      first <- columns[[p]][1]
      last <- first - 1 + length(smaller$columns[[p]]) 
      for (i in 1:length(smaller)) {
        new.smaller.m[i,first:last] <- smaller$m[i,smaller$columns[[p]]]
      }
    }
  } else {
    smaller <- c1
    bigger <- c2
    new.smaller.m <- smaller$m
  }
  bigger$m <- rbind(new.smaller.m, bigger$m)
  return(bigger)
}

average.qchain <- function(ch) {
  params <- names(ch$columns)
  len <- 0
  for (p in params) {
    len <- len + length(ch$columns[[p]])
  }
  m <- rep(0, len)
  for (p in params) {
    i <- ch$columns[[p]]
    m[i] <- average(attr(ch$model.fn, "types")[[p]](ch$m[,i,drop=F]))
  }
  return(ch$model.fn(ch$data, ch$columns, m, off.by.one=ch$off.by.one))
}

DIC.qchain <- function(ch) {
  # Compute deviance for point estimator
  dex <- bdev(average(ch))
  # Average deviance for models in this chain
  exd <- 0
  for (i in 1:length(ch)) {
    exd <- exd + bdev(ch[[i]])
  }
  exd <- exd/length(ch)
  # Return DIC
  o <- list(dex=dex, exd=exd, DIC=(2*exd-dex))
  return(o)
}

effective.params.qchain <- function(ch) {
  dex <- bdev(average(ch))
  return((DIC(ch)-dex)/2)
}

save.qchain <- function(ch, name) {
  datestr <- format(Sys.time(), "%Y_%m_%d")
  filename <- paste(paste(name, datestr, sep='_'), ".RData", sep='')
  remove.old <- F
  if (file.exists(filename)) { 
    old.filename <- paste(filename, ".old", sep='')
    file.rename(filename, old.filename)
    remove.old <- T
  }
  base:::save(ch, file=filename)
  if (remove.old) {
    if (file.exists(old.filename)) {
      file.remove(old.filename)
    }
  }
}

mixture.qchain <- function(ch) {
  m <- average(ch)
  m$raw <- ch
#  m$stats <- list()
#    DIC=function() {
#      DIC(m$raw)$DIC
#    },
#    p=function() {
#      bayesianp(m)
#    K=function() {
#      length(m$active)
#    }
#  )
#  environment(m$stats$K) <- emptyenv()
  return(m)
}

mixture.jobjRef <- function(chobj, data, model.fn, jqgmodel) {
  ch <- qchain(chobj, data, model.fn, jqgmodel)
  if (length(.jmethods(J(jqgmodel), "pointEstimate")) > 0) {
    ptechlobj <- J(jqgmodel)$pointEstimate(chobj, jqgobj(data$data))
    ptechobj <- .jnew("org/jqgibbs/Chain")
    ptechobj$addLink(ptechlobj)
    ptech <- qchain(ptechobj, data, model.fn, jqgmodel)
    pte <- ptech$m[1,]
    m <- model.fn(data, ptech$columns, pte, off.by.one=ptech$off.by.one)
  } else {
    m <- average(ch)
  }
  m$raw <- ch
  return(m)
}

acor.qchain <- function(ch, tau=1) {
  m <- ch$m
  nccol <- ncol(m)
  mn <- colMeans(m)
  n <- t(t(m) - mn)
  acmat <- matrix(0, nrow(m)-tau, nccol)
  varmat <- matrix(0, nrow(m)-tau, nccol)
  acmat[1,] <- n[1,]*n[1+tau,]
  varmat[1,] <- n[1,]^2
  if (nrow(m) > 2) {
    for (j in 2:(nrow(m)-tau)) {
      acmat[j,] <- n[j,]*n[j+tau,]
      varmat[j,] <- n[j,]^2
    }
  }
  ac <- colSums(acmat)
  vr <- colSums(varmat)
  vr[vr==0] <- NA
  return(ac/vr)
}

###########################################################################
# jqgibbs Experiment subclass
###########################################################################

jqgexpt <- function(name, setup) {
  e <- expt(name, setup)
  e$run <- function(data, obs) {
    environment(e$setup) <- new.env()
    settings <- e$setup(data, obs.to.predictors(obs))
    m <- jqgibbs(data, settings$model.fn, settings$jqgmodel, name, settings$samplerClass,
            settings$hypers, settings$init, settings$length, settings$burnin,
            settings$lag, settings$savetime, settings$finalsave)
    return(m)
  }
  class(e) <- c(class(e), "jqgexpt")
  return(e)
}

