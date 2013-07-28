##############################################################################
# dlp.R
# Functions for running DLP experiments.
##############################################################################
# 16/2/2011: Created
# 20/2/2011: Added functions for splitting up results files
##############################################################################
# Copyright Ewan Dunbar
# University of Maryland
##############################################################################

##############################################################################
# Generic methods
##############################################################################
resultid <- function(...) {
  UseMethod("resultid")
}

##############################################################################
# Load a DLP dataset
##############################################################################
dataset.character <- function(dsname) {
  if (dsname %in% data(package="dlp")$results[,"Item"]) {
    data(list=dsname, package="dlp")
    d <- dataset(get(dsname))
    rm(dsname)
  } else if (dsname %in% ls(".GlobalEnv")) {
    warning("Loading dataset from global environment")
    d <- get(dsname)
  }
  return(d)
}

##############################################################################
# Write a DLP dataset.
#
# Dataset is located in clustering.R, but this method is overloaded
# here because we have a somewhat idiosyncratic way of keeping track
# of data (txt+R). This might change in the future.
#
# This function really has no impact currently,
# since the only purpose of the txt+R format is to change the column
# names so that the data frame can be converted to a dataset, but the
# output of as.data.frame.dataset has column names which are
# exactly what they should be already. The function should be used, however,
# in case this state of affairs changes.
##############################################################################
write.dataset <- function(ds, name) {
  # Write data frame
  dfds <- as.data.frame.dataset(ds)
  dfcolnames <- default.colnames(ds)
  fndf <- paste(name, ".txt", sep="")
  write.table(dfds, file=fndf, row.names=F, col.names=T)
  # Write R code to load the data frame
  rcode <- c(
    paste(name, " <- read.table(\"", fndf, "\", header=T)", sep=''),
    paste("names(", name, ") <- c(\"", paste(dfcolnames, collapse='\",\"'), "\")", sep='')
  )
  fnrc <- paste(name, ".R", sep="")
  cat(rcode, file=fnrc, sep='\n')
}

##############################################################################
# A handy function for generating the "corrected" data sets for
# the "classical" DLP experiments (manual correction). Will only be
# correct for one predictor currently. All other predictors will remain,
# but the one being corrected for will be removed.
##############################################################################
correct <- function(dat, byC=F) {
  ds <- dataset(dat)
  X <- ds$data
  O <- ds$tronly
  C <- ds$classes

  corrds <- ds
  if (byC) {
    for (cl in levels(C)) {
      xin <- X[C == cl & O[,1] == 1,]
      xelse <- X[C == cl & O[,1] == 0,]
      correction <- colMeans(xin) - colMeans(xelse)
      for (i in 1:ncol(X)) {
        corrds$data[C == cl & O[,1] == 1,i] <- xin[,i] - correction[i]
      }
    }
  } else {
    xin <- X[O[,1] == 1,]
    xelse <- X[O[,1] == 0,]
    correction <- colMeans(xin) - colMeans(xelse)
    for (i in 1:ncol(X)) {
      corrds$data[O[,1] == 1,i] <- xin[,i] - correction[i]
    }
  }
  corrds$tronly <- O[,-1,drop=F]
  return(corrds)
}

##############################################################################
# Cross-validation
##############################################################################
folds <- function(C, k, maintainMix=T) {
  f <- as.list(rep(NA, k))
  if (maintainMix) {
    classes <- names(summary(as.factor(C)))
    for (cl in classes) {
      ci <- (1:length(C))[C==cl]
      fci <- folds(ci, k, maintainMix=F)
      for (i in 1:k) {
        fck <- ci[fci[[i]]]
        if (!is.na(f[[i]][1])) {
          f[[i]] <- c(f[[i]], fck)
        } else {
          f[[i]] <- fck
        }
      }
    }
  } else {
    si <- sample(1:length(C))
    endi <- 0
    nCdk <- floor(length(C)/k)
    for (i in 1:k) {
      starti <- endi + 1
      endi <- starti + nCdk-1
      f[[i]] <- si[starti:endi]
    }
  }
  return(f)
}

kfxv <- function(expt, ds, k, maintainMix, savefile=NULL, fold.indices=NULL) {
  X <- ds$data
  O <- center.obs(ds$tronly)
  C <- ds$classes
  if (identical(fold.indices, NULL)) {
    Oc <- apply(O, 1, function(x) paste(x, collapse=''))
    CO <- paste(as.character(C), Oc, sep='')
    fold.indices <- folds(CO, k, maintainMix)
  }
  r <- list()
  r$byfold <- list()
  for (i in 1:length(fold.indices)) {
    # Get training columns, test columns, and class vector
    # for this fold
    test.indices <- fold.indices[[i]]
    trng.indices <- (1:nrow(X))[!(1:nrow(X) %in% fold.indices[[i]])]
    XTi <- X[trng.indices,,drop=F]
    Oi <- O[trng.indices,,drop=F]
    Xi <- X[test.indices,,drop=F]
    Ci <- C[test.indices]
    # Perform clustering
    write(paste("Running fold", i, "\n"), file="")
    model <- expt$run(XTi, Oi)
    # Save results
    r$byfold[[i]] <- list()
    r$byfold[[i]]$model <- model
    r$byfold[[i]]$ds <- ds
    r$byfold[[i]]$trng.indices <- trng.indices
    # Write to savefile
    if (!identical(savefile, NULL)) {
      save(r, file=savefile)
    }
  }
  # Return
  return(r)
}

########################################################################
# Convenient function for running experiments. Not that clever, as it only
# saves after all the cross-validation for a single experiment is done, and
# even then it doesn't save incrementally or free any memory.
########################################################################
run_experiments <- function(pairs, filename, kfxv_save_file=NULL, fold.indices=NULL) {
  results <- list()
  index <- c()
  i <- 1
  for (pair in pairs) {
    # Load dataset
    ds <- pair[1]
    write(paste("Data set '", ds, "'\n", sep=""), file="")
    d <- dataset(ds)
    # Run experiments
    expt <- experiments[[pair[2]]]
    write(paste("Experiment '", expt$name, "'\n", sep=""), file="")
    # Run experiment
    r <- kfxv(expt, d, k=10, maintainMix=T, savefile=kfxv_save_file, fold.indices)
    # Save results
    index[i] <- paste(pair[1], pair[2], sep="_")
    results[[i]] <- list(expt, r)
    i <- i+1
    print("Saving image...")
    save(index, results, file=filename)
  }
}

########################################################################
# Structure holding an identifier for a set of results
########################################################################
resultid.character <- function(s, id_table=NULL, level="foldrun") {
  # Create empty resultid
  rid <- list()
  # Split string
  fields <- strsplit(s, "_")[[1]]
  # Check to see if there is an initial X field - this was used ad hoc
  # to distinguish separate exptruns; here we just strip it and let
  # the check against id_table do its work
  if ((substr(fields[1], 1, 1) == "X") & 
      (nchar(fields[1]) > 1) &
      (grepl("^[[:digit:]]+$", substr(fields[1], 2, nchar(fields[1]))))) {
      fields <- fields[-1]
  }
  # Find the experiment name
  en_fields <- fields %in% names(experiments)
  if (sum(en_fields) != 1) {
    rid$experiment <- NULL
  } else {
    rid$experiment <- fields[en_fields][1]  
    fields <- fields[-(which(en_fields)[1])]
  }
  # Fill in run info, if present
  rid$run <- NULL
  rid$fold <- NULL
  rid$foldrun <- NULL
  remove <- c()
  for (i in 1:length(fields)) {
    f <- fields[i]
    if (grepl("^[[:digit:]]+$", f)) {
      if (identical(rid$run, NULL)) {
        rid$run <- as.numeric(f)
        remove <- c(remove, i)
      } else if (identical(rid$fold, NULL)) {
        rid$fold <- as.numeric(f)
        remove <- c(remove, i)
      } else if (identical(rid$foldrun, NULL)) {
        rid$foldrun <- as.numeric(f)
        remove <- c(remove, i)
      }
    }
  }
  if (!(identical(remove, NULL))) {
    fields <- fields[-1*remove]
  }
  # The rest should be a dataset
  ds <- paste(fields, collapse="_")
  if (ds %in% data(package="dlp")$results[,"Item"]) {
    rid$ds <- ds
  } else if (ds %in% ls(".GlobalEnv")) {
    rid$ds <- ds
    warning("Assuming dataset from global environment")
  } else {
    rid$ds <- NULL
  }
  class(rid) <- "resultid"
  # Fill in remaining fields from id_table
  if (identical(rid$run, NULL) |
      identical(rid$fold, NULL) |
      identical(rid$foldrun, NULL)) {
    rid <- resultid(rid, id_table, level)
  }
  # Return resultid object
  return(rid)
}

# FIXME - What a monster!
resultid.resultid <- function(rid, id_table=NULL, level="foldrun") {
  if (level == "run") {
    last_run <- 0
    if (!(identical(id_table, NULL))) {
      previous_entries <- id_table[id_table$EXPT == rid$experiment &
                                   id_table$DATA == rid$ds,]
      if (nrow(previous_entries) > 0) {
        last_run <- max(previous_entries$RUN)
      }
    }
    rid$run <- last_run + 1
    rid$fold <- NULL
    rid$foldrun <- NULL
    return(rid)
  }
  if (level == "fold") {
    last_run <- 1
    last_fold <- 0
    if (!(identical(id_table, NULL))) {
      previous_runs <- id_table[id_table$EXPT == rid$experiment &
                                id_table$DATA == rid$ds,]
      if (identical(rid$run, NULL)) {
        last_run <- max(previous_runs$RUN)
        previous_entries <- previous_runs[previous_runs$RUN == last_run,]
      } else {
        previous_entries <- previous_runs[previous_runs$RUN == rid$run,]
      }
      if (nrow(previous_entries) > 0) {
        last_fold <- max(previous_entries$FOLD)
      }
    }
    if (identical(rid$run, NULL)) {
      rid$run <- last_run
    }
    rid$fold <- last_fold + 1
    rid$foldrun <- NULL
    return(rid)
  }
  if (level == "foldrun") {
    last_run <- 1
    last_fold <- 1
    last_foldrun <- 0
    if (!(identical(id_table, NULL))) {
      previous_runs <- id_table[id_table$EXPT == rid$experiment &
                                id_table$DATA == rid$ds,]
      if (identical(rid$run, NULL)) {
        last_run <- max(previous_runs$RUN)
        previous_folds <- previous_runs[previous_runs$RUN == last_run,]
      } else {
        previous_folds <- previous_runs[previous_runs$RUN == rid$run,]
      }
      if (identical(rid$fold, NULL)) {
        last_fold <- max(previous_folds$FOLD)
        previous_entries <- previous_folds[previous_folds$FOLD == last_fold,]
      } else {
        previous_entries <- previous_folds[previous_folds$FOLD == rid$fold,]
      }
      if (nrow(previous_entries) > 0) {
        last_foldrun <- max(previous_entries$FOLDRUN)
      }
    }
    if (identical(rid$run, NULL)) {
      rid$run <- last_run
    }
    if (identical(rid$fold, NULL)) {
      rid$fold <- last_fold
    }
    rid$foldrun <- last_foldrun + 1
    return(rid)
  }
  return(NULL)
}

as.character.resultid <- function(id) {
  ns <- c()
  if (!identical(id$ds, NULL)) {
    ns <- c(ns, id$ds)
  }
  if (!identical(id$experiment, NULL)) {
    ns <- c(ns, id$experiment)
  }
  if (!identical(id$run, NULL)) {
    ns <- c(ns, as.character(id$run))
  }
  if (!identical(id$fold, NULL)) {
    ns <- c(ns, as.character(id$fold))
  }
  if (!identical(id$foldrun, NULL)) {
    ns <- c(ns, as.character(id$foldrun))
  }
  if (!identical(ns, NULL)) {
    return(paste(ns, collapse="_"))
  } else {
    stop("Cannot generate string from empty id")
  }
}

complete.resultid <- function(rid) {
  if (identical(rid$experiment, NULL)) {
    return(FALSE)
  }
  if (identical(rid$ds, NULL)) {
    return(FALSE)
  }
  if (identical(rid$run, NULL)) {
    return(FALSE)
  }
  if (identical(rid$fold, NULL)) {
    return(FALSE)
  }
  if (identical(rid$foldrun, NULL)) {
    return(FALSE)
  }
  return(TRUE)
}

########################################################################
# The files stored by run_experiments and kfxv are somewhat unwieldy.
# Rather than rewriting these functions now, I have provided the following
# function to store reduced versions of the models in these files. It
# also backs up the results that were stored, although see the following
# function for a more convenient way to recompute statistics from stored
# model files.
########################################################################
upgrade_savefile <- function(savefile, id_table, max_folds=10, dir=".", use_rownames=F) {
  # Load the savefile
  load(savefile)
  # Check to make sure the savefile has an index and results
  if (!(exists("index", inherits=F) &
        exists("results", inherits=F) &
       (length(index) == length(results)))) {
    stop("Missing or inconsistent index or results")
  }
  # For each entry in the index
  i_results <- 1
  for (label in index) {
    # Extract the id from the label
    id <- resultid(label, id_table, level="run")
    if (identical(id$experiment, NULL) | identical(id$ds, NULL)) {
      warning(paste("Unknown experiment or dataset in", label))
    } else {
      # Separate experiment details, folds, and summary
      expt <- results[[i_results]][[1]]
      folds <- results[[i_results]][[2]][[1]]
      max_sr_index <- length(results[[i_results]][[2]])
      summ <- results[[i_results]][[2]][2:max_sr_index]
      # Back up the previously obtained summary
      sr_fn <- paste(as.character(id), ".old_summary.RData", sep='')
      sr_fn <- paste(dir, sr_fn, sep='/')
      save(summ, expt, file=sr_fn)
      # For each element in the list
      for (fold in results[[i_results]][[2]][[1]]) {
        # Fill in the complete id
        id <- resultid(id, id_table, level="fold")
        id$foldrun <- 1
        if ((identical(id_table, NULL) || !(as.character(id) %in% id_table$LABEL)) &&
            (id$fold <= max_folds)) {
          # Deduce the appropriate training and test sets
          if (use_rownames) {
            orig <- dataset(id$ds)
            fulln <- rownames(orig$data)
            trn <- rownames(fold$model$data)
            testn <- fulln[!(fulln %in% trn)]
            dss <- list()
            dss$training <- dataset(as.data.frame(dataset(id$ds))[trn,,drop=F])
            dss$test <- dataset(as.data.frame(dataset(id$ds))[testn,,drop=F])
          } else {
            dss <- find_training_test(fold$model$data, id$ds)
          }
          save(dss, file="SAVED_TRTEST.RData")
          if (identical(dss, NULL)) {
            warning(paste("Unable to find training/test data for", id))
          } else {
            m <- fold$model
            # Eliminate junk
            if (length(m$stats) > 0) {
              for (i in 1:length(m$stats)) {
                environment(m$stats[[i]]) <- emptyenv()
              }
            }
            # Save the full raw model
            raw_fn <- paste(as.character(id), ".raw.RData", sep='')
            raw_fn <- paste(dir, raw_fn, sep='/')
            save(m, dss, expt, file=raw_fn)
            # Save the point estimate model
            m$raw <- NULL
            pte_fn <- paste(as.character(id), ".pte.RData", sep='')
            pte_fn <- paste(dir, pte_fn, sep='/')
            save(m, dss, expt, file=pte_fn)
            # Update the id table
            id_table <- update_idtable(id, id_table)
          }
        }
      }
    }
    i_results <- i_results + 1
  }
  rm(results)
  rm(index)
  gc()
  return(id_table)
}

find_training_test <- function(tr, dsname) {
  orig_data <- dataset(dsname)$data
  if (nrow(tr) > nrow(orig_data)) {
    warning("Data set is too small to be training data source")
    return(NULL)
  }
  if ((ncol(tr) != ncol(orig_data)) |
      !identical(colnames(tr), colnames(orig_data))) {
    warning("Data set has mismatching columns")
    if (prod(colnames(tr) %in% colnames(orig_data)) == 1) {
      orig_data <- orig_data[,colnames(tr),drop=F]
    } else {
      return(NULL)
    }
  }
  tr_m <- as.matrix(tr)
  orig_data_m <- as.matrix(orig_data)
  
  orig_by_tr_r <- sort(rep(1:nrow(orig_data_m), nrow(tr_m)))
  orig_by_tr <- orig_data_m[orig_by_tr_r,,drop=F]
  tr_rep_r <- rep(1:nrow(tr_m), nrow(orig_data_m))
  tr_rep <- tr_m[tr_rep_r,,drop=F]

  candidates <- rep(T, nrow(orig_by_tr))
  equal <- matrix(0, nrow(orig_by_tr), ncol(orig_by_tr))
  for (i in 1:ncol(orig_by_tr)) {
    equal[candidates,i] <- orig_by_tr[candidates,i] == tr_rep[candidates,i]
    if (i < ncol(orig_by_tr)) {
      candidates <- as.logical(equal[,i])
    }
  }
  equal_all <- equal[,ncol(orig_by_tr)]
  equal_mat <- matrix(equal_all, nrow(tr_m), nrow(orig_data_m))

  while (prod(rowSums(equal_mat) > 1) == 1) {
    dupe_row_num <- which(rowSums(equal_mat) > 1)[1]
    dupe_row <- equal_mat[dupe_row_num,]
    dupe_col_nums <- which(as.logical(dupe_row))
    dupe_col_sums <- colSums(equal_mat[-dupe_row_num,dupe_col_nums])
    if (sum(dupe_col_sums>0) >= 1) {
      first_available <- which(dupe_col_sums>0)[1]
    } else {
      first_available <- dupe_col_nums[1]
    }
    dupe_row[first_available] <- 0
  }

  orig_tr_indices <- which(as.logical(colSums(equal_mat)))

  if (length(orig_tr_indices) != nrow(tr)) {
    warning("Data set does not match this training data")
    return(NULL)
  }

  dss <- list()
  dss$training <- dataset(as.data.frame(dataset(dsname))[orig_tr_indices,,drop=F])
  dss$test <- dataset(as.data.frame(dataset(dsname))[!(1:nrow(orig_data) %in% orig_tr_indices),,drop=F])
  return(dss)
}

update_idtable <- function(id, id_table) {
  id_row <- data.frame(LABEL=as.character(id),
                         EXPT=id$experiment,
                         DATA=id$ds,
                         RUN=id$run,
                         FOLD=id$fold,
                         FOLDRUN=id$foldrun)
  if (identical(id_table, NULL)) {
    return(id_row)
  } else {
    return(rbind(id_table, id_row))
  }
}

########################################################################
# Run statistics on a result
########################################################################
run_stats <- function(resultid, result_table, stats, dir='.') {
  load(paste(dir, "/", as.character(resultid), ".pte.RData", sep=''))
  res_row <- data.frame(LABEL=as.character(resultid))
  for (f in stats) {
    res_row <- cbind(res_row, f(m, dss))
  }
  if (!identical(names(stats), NULL)) {
    names(res_row)[2:ncol(res_row)] <- names(stats)
  }
  if (identical(result_table, NULL)) {
    return(res_row)
  } else {
    return(rbind(result_table, res_row))
  }
}

pairwise_prec.dlp.tr <- function(m, dss) {
  return(precision(m$Z, dss$training$classes, pairwise=T))
}

pairwise_prec.dlp.te <- function(m, dss) {
  prediction <- classify(m, dss$test$data)
  return(precision(prediction, dss$test$classes, pairwise=T))
}

pairwise_prec.prod.dlp.tr <- function(m, dss) {
  if (ncol(dss$training$tronly) > 0) {
    truth_ch <- as.character(dss$training$classes)
    for (i in 1:ncol(dss$training$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$training$tronly[,i]))
    }
    return(precision(m$Z, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_prec.prod.dlp.te <- function(m, dss) {
  if (ncol(dss$test$tronly) > 0) {
    prediction <- classify(m, dss$test$data)
    truth_ch <- as.character(dss$test$classes)
    for (i in 1:ncol(dss$test$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$test$tronly[,i]))
    }
    return(precision(prediction, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_rec.dlp.tr <- function(m, dss) {
  return(recall(m$Z, dss$training$classes, pairwise=T))
}

pairwise_rec.dlp.te <- function(m, dss) {
  prediction <- classify(m, dss$test$data)
  return(recall(prediction, dss$test$classes, pairwise=T))
}

pairwise_rec.prod.dlp.tr <- function(m, dss) {
  if (ncol(dss$training$tronly) > 0) {
    truth_ch <- as.character(dss$training$classes)
    for (i in 1:ncol(dss$training$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$training$tronly[,i]))
    }
    return(recall(m$Z, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_rec.prod.dlp.te <- function(m, dss) {
  if (ncol(dss$test$tronly) > 0) {
    prediction <- classify(m, dss$test$data)
    truth_ch <- as.character(dss$test$classes)
    for (i in 1:ncol(dss$test$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$test$tronly[,i]))
    }
    return(recall(prediction, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_f.dlp.tr <- function(m, dss) {
  return(fscore(m$Z, dss$training$classes, pairwise=T))
}

pairwise_f.dlp.te <- function(m, dss) {
  prediction <- classify(m, dss$test$data)
  return(fscore(prediction, dss$test$classes, pairwise=T))
}

pairwise_f.prod.dlp.tr <- function(m, dss) {
  if (ncol(dss$training$tronly) > 0) {
    truth_ch <- as.character(dss$training$classes)
    for (i in 1:ncol(dss$training$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$training$tronly[,i]))
    }
    return(fscore(m$Z, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}


pairwise_f.prod.dlp.te <- function(m, dss) {
  if (ncol(dss$test$tronly) > 0) {
    prediction <- classify(m, dss$test$data)
    truth_ch <- as.character(dss$test$classes)
    for (i in 1:ncol(dss$test$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$test$tronly[,i]))
    }
    return(fscore(prediction, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_f.alpr.dlp.tr <- function(m, dss) {
  if ((ncol(dss$training$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    prediction <- classify(m, dss$training$data, predict_B=T)
    truth_ch <- as.character(dss$training$classes)
    for (i in 1:ncol(dss$training$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$training$tronly[,i]))
    }
    return(fscore(prediction, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_f.alpr.dlp.te <- function(m, dss) {
  if ((ncol(dss$training$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    prediction <- classify(m, dss$test$data, predict_B=T)
    truth_ch <- as.character(dss$test$classes)
    for (i in 1:ncol(dss$test$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$test$tronly[,i]))
    }
    return(fscore(prediction, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_prec.alpr.dlp.tr <- function(m, dss) {
  if ((ncol(dss$training$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    prediction <- classify(m, dss$training$data, predict_B=T)
    truth_ch <- as.character(dss$training$classes)
    for (i in 1:ncol(dss$training$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$training$tronly[,i]))
    }
    return(precision(prediction, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_prec.alpr.dlp.te <- function(m, dss) {
  if ((ncol(dss$training$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    prediction <- classify(m, dss$test$data, predict_B=T)
    truth_ch <- as.character(dss$test$classes)
    for (i in 1:ncol(dss$test$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$test$tronly[,i]))
    }
    return(precision(prediction, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_rec.alpr.dlp.tr <- function(m, dss) {
  if ((ncol(dss$training$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    prediction <- classify(m, dss$training$data, predict_B=T)
    truth_ch <- as.character(dss$training$classes)
    for (i in 1:ncol(dss$training$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$training$tronly[,i]))
    }
    return(recall(prediction, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_rec.alpr.dlp.te <- function(m, dss) {
  if ((ncol(dss$training$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    prediction <- classify(m, dss$test$data, predict_B=T)
    truth_ch <- as.character(dss$test$classes)
    for (i in 1:ncol(dss$test$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$test$tronly[,i]))
    }
    return(recall(prediction, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}

K.dlp.tr <- function(m, dss) {
  return(length(unique(m$Z)))
}

K.dlp.te <- function(m, dss) {
  prediction <- classify(m, dss$test$data)
  return(length(unique(prediction)))
}

pairwise_prec.dlp.te.fb <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B)
    return(precision(prediction, dss$test$classes, pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_prec.prod.dlp.te.fb <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B)
    truth_ch <- as.character(dss$test$classes)
    for (i in 1:ncol(dss$test$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$test$tronly[,i]))
    }
    return(precision(prediction, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_rec.dlp.te.fb <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B)
    return(recall(prediction, dss$test$classes, pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_rec.prod.dlp.te.fb <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B)
    truth_ch <- as.character(dss$test$classes)
    for (i in 1:ncol(dss$test$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$test$tronly[,i]))
    }
    return(recall(prediction, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_f.dlp.te.fb <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B)
    return(fscore(prediction, dss$test$classes, pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_f.prod.dlp.te.fb <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B)
    truth_ch <- as.character(dss$test$classes)
    for (i in 1:ncol(dss$test$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$test$tronly[,i]))
    }
    return(fscore(prediction, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_f.alpr.dlp.te.fb <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, predict_B=T, fixed_B=B)
    truth_ch <- as.character(dss$test$classes)
    for (i in 1:ncol(dss$test$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$test$tronly[,i]))
    }
    return(fscore(prediction, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_prec.alpr.dlp.te.fb <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, predict_B=T, fixed_B=B)
    truth_ch <- as.character(dss$test$classes)
    for (i in 1:ncol(dss$test$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$test$tronly[,i]))
    }
    return(precision(prediction, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}

pairwise_rec.alpr.dlp.te.fb <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, predict_B=T, fixed_B=B)
    truth_ch <- as.character(dss$test$classes)
    for (i in 1:ncol(dss$test$tronly)) {
      truth_ch <- paste(truth_ch, as.character(dss$test$tronly[,i]))
    }
    return(recall(prediction, factor(truth_ch), pairwise=T))
  } else {
    return(NA)
  }
}

K.dlp.te.fb <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B)
    return(length(unique(prediction)))
  } else {
    return(NA)
  }
}

pairwise_prec.sec1.dlp.tr <- function(m, dss) {
  return(precision(m$Z, dss$training$secclasses[,1], pairwise=T))
}

pairwise_prec.sec1.dlp.te.fbif <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B)
  } else {
    prediction <- classify(m, dss$test$data)
  }
  return(precision(prediction, dss$test$secclasses[,1], pairwise=T))
}

pairwise_prec.sec2.dlp.tr <- function(m, dss) {
  return(precision(m$Z, dss$training$secclasses[,2], pairwise=T))
}

pairwise_prec.sec2.dlp.te.fbif <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B)
  } else {
    prediction <- classify(m, dss$test$data)
  }
  return(precision(prediction, dss$test$secclasses[,2], pairwise=T))
}

pairwise_prec.dlp.te.fbif <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B)
  } else {
    prediction <- classify(m, dss$test$data)
  }
  return(precision(prediction, dss$test$classes, pairwise=T))
}


pairwise_prec.dlp.te4.fbif <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B, max=4)
  } else {
    prediction <- classify(m, dss$test$data, max=4)
  }
  return(precision(prediction, dss$test$classes, pairwise=T))
}

pairwise_rec.sec1.dlp.tr <- function(m, dss) {
  return(recall(m$Z, dss$training$secclasses[,1], pairwise=T))
}

pairwise_rec.sec1.dlp.te.fbif <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B)
  } else {
    prediction <- classify(m, dss$test$data)
  }
  return(recall(prediction, dss$test$secclasses[,1], pairwise=T))
}

pairwise_rec.sec2.dlp.tr <- function(m, dss) {
  return(recall(m$Z, dss$training$secclasses[,2], pairwise=T))
}

pairwise_rec.sec2.dlp.te.fbif <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B)
  } else {
    prediction <- classify(m, dss$test$data)
  }
  return(recall(prediction, dss$test$secclasses[,2], pairwise=T))
}

pairwise_rec.dlp.te.fbif <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B)
  } else {
    prediction <- classify(m, dss$test$data)
  }
  return(recall(prediction, dss$test$classes, pairwise=T))
}

pairwise_rec.dlp.te4.fbif <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B, max=4)
  } else {
    prediction <- classify(m, dss$test$data, max=4)
  }
  return(recall(prediction, dss$test$classes, pairwise=T))
}

pairwise_f.sec1.dlp.tr <- function(m, dss) {
  return(fscore(m$Z, dss$training$secclasses[,1], pairwise=T))
}

pairwise_f.sec1.dlp.te.fbif <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B)
  } else {
    prediction <- classify(m, dss$test$data)
  }
  return(fscore(prediction, dss$test$secclasses[,1], pairwise=T))
}

pairwise_f.sec2.dlp.tr <- function(m, dss) {
  return(fscore(m$Z, dss$training$secclasses[,2], pairwise=T))
}

pairwise_f.sec2.dlp.te.fbif <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B)
  } else {
    prediction <- classify(m, dss$test$data)
  }
  return(fscore(prediction, dss$test$secclasses[,2], pairwise=T))
}

pairwise_f.dlp.te.fbif <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B)
  } else {
    prediction <- classify(m, dss$test$data)
  }
  return(fscore(prediction, dss$test$classes, pairwise=T))
}

pairwise_f.dlp.te4.fbif <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B, max=4)
  } else {
    prediction <- classify(m, dss$test$data, max=4)
  }
  return(fscore(prediction, dss$test$classes, pairwise=T))
}


K.dlp.te.fbif <- function(m, dss) {
  if ((ncol(dss$test$tronly) > 0) & ("flgfa" %in% class(m)) && (ncol(m$B) > 1)) {
    B <- cbind(1, dss$test$tronly)
    prediction <- classify(m, dss$test$data, fixed_B=B)
    return(length(unique(prediction)))
  } else {
    return(K.dlp.te(m, dss))
  }
}

sec1.stats <- function() {
  l <- list(K.dlp.te.fbif, K.dlp.tr, 
                                     pairwise_f.sec1.dlp.te.fbif,
                                     pairwise_f.sec1.dlp.tr,
                                     pairwise_prec.sec1.dlp.te.fbif,
                                     pairwise_prec.sec1.dlp.tr,
                                     pairwise_rec.sec1.dlp.te.fbif,
                                     pairwise_rec.sec1.dlp.tr,
                                     pairwise_f.dlp.te.fbif,
                                     pairwise_f.dlp.tr,
                                     pairwise_prec.dlp.te.fbif,
                                     pairwise_prec.dlp.tr,
                                     pairwise_rec.dlp.te.fbif,
                                     pairwise_rec.dlp.tr)

  names(l) <- c("Ktefbi", "Ktr",
                                "pFte1fbi", "pF1tr",
                                "pprte1fbi", "ppr1tr",
                                "prete1fbi", "pre1tr",
                                "pFtefbi", "pFtr",
                                "pprtefbi", "pprtr",
                                "pretefbi", "pretr")
  return(l)

}

sec2.stats <- function() {
  l <- list(
                                     pairwise_f.sec2.dlp.te.fbif,
                                     pairwise_f.sec2.dlp.tr,
                                     pairwise_prec.sec2.dlp.te.fbif,
                                     pairwise_prec.sec2.dlp.tr,
                                     pairwise_rec.sec2.dlp.te.fbif,
                                     pairwise_rec.sec2.dlp.tr)

  names(l) <- c(
                                "pFte2fbi", "pF2tr",
                                "pprte2fbi", "ppr2tr",
                                "prete2fbi", "pre2tr")
  return(l)

}

standard.stats <- function() {
  l <- list(K.dlp.te, K.dlp.tr,
            pairwise_f.dlp.te, pairwise_prec.dlp.te, pairwise_rec.dlp.te,
            pairwise_f.dlp.tr, pairwise_prec.dlp.tr, pairwise_rec.dlp.tr,
            pairwise_f.prod.dlp.te, pairwise_prec.prod.dlp.te, pairwise_rec.prod.dlp.te,
            pairwise_f.prod.dlp.tr, pairwise_prec.prod.dlp.tr, pairwise_rec.prod.dlp.tr,
            pairwise_f.alpr.dlp.te, pairwise_prec.alpr.dlp.te, pairwise_rec.alpr.dlp.te,
            pairwise_f.alpr.dlp.tr, pairwise_prec.alpr.dlp.tr, pairwise_rec.alpr.dlp.tr)
  names(l) <- c("Kte", "Ktr",
                              "pFte", "pprte", "prete",
                              "pFtr", "pprtr", "pretr",
                              "pFprodte", "pprprodte", "preprodte",
                              "pFprodtr", "pprprodtr", "preprodtr",
                              "pFalprte", "ppralprte", "prealprte",
                              "pFalprtr", "ppralprtr", "prealprtr")
  return(l)
}

fb.stats <- function() {
  l <- list(K.dlp.te.fb,
            pairwise_f.dlp.te.fb, pairwise_prec.dlp.te.fb, pairwise_rec.dlp.te.fb,
            pairwise_f.prod.dlp.te.fb, pairwise_prec.prod.dlp.te.fb, pairwise_rec.prod.dlp.te.fb,
            pairwise_f.alpr.dlp.te.fb, pairwise_prec.alpr.dlp.te.fb, pairwise_rec.alpr.dlp.te.fb)
  names(l) <- c("Ktefb", 
                              "pFtefb", "pprtefb", "pretefb",
                              "pFprodtefb", "pprprodtefb", "preprodtefb",
                              "pFalprtefb", "ppralprtefb", "prealprtefb")
  return(l)

}

fbi.stats <- function() {
  l <- list(K.dlp.tr, K.dlp.te.fbif,
            pairwise_f.dlp.tr, pairwise_prec.dlp.tr, pairwise_rec.dlp.tr,
            pairwise_f.dlp.te.fbif, pairwise_prec.dlp.te.fbif, pairwise_rec.dlp.te.fbif,
            pairwise_f.dlp.te4.fbif, pairwise_prec.dlp.te4.fbif, pairwise_rec.dlp.te4.fbif)
  names(l) <- c("Ktr", "Ktefbi", 
                              "pFtr", "pprtr", "pretr",
                              "pFtefbi", "pprtefbi", "pretefbi",
                              "pFte4fbi", "pprte4fbi", "prete4fbi")
  return(l)

}
###########################################################################
# Special plotting functions
###########################################################################

#dlp.transform <- function(d) {
#  d[,2:1]
#}
#
#dlp.initplot <- function(m, pretty=T) {
#  if (pretty) {
#    par(bg="#fcffdd", family="Georgia")
#  }
#  K <- length(unique(m$Z))
#  r <- rainbow(3*K)
#  r[1:K] <- r[3*(0:(K-1))+1]
#  palette(r)
#  palette(adjustcolor(palette(), alpha.f=0.25))
#  plot.new()
#  xlim <- default.plotlim(dlp.transform(m$data$data)[,1])[2:1]
#  ylim <- default.plotlim(dlp.transform(m$data$data)[,2])[2:1]
#  plot.window(xlim, ylim)
#  par(xaxp=c(xlim,5), yaxp=c(ylim,5))
#  axis(1)
#  axis(2)
#}
#
#dlp.plot <- function(m, ...) {
#  plot(m, initplot=dlp.initplot, transform=dlp.transform, pch=20, ...)
#}
#
#paper.initplot <- function (m) {
#    par(bg = "#ffffff", family = "Times")
#    K <- length(unique(m$Z))
#    r <- rainbow(3*K, v=0)
#    r[1:K] <- r[3*(0:(K-1))+1]
#    palette(r)
#    palette(adjustcolor(palette(), alpha.f = 0.25))
#    plot.new()
#    xlim <- default.plotlim(dlp.transform(m$data)[, 1])[2:1]
#    ylim <- default.plotlim(dlp.transform(m$data)[, 2])[2:1]
#    plot.window(xlim, ylim)
#    par(xaxp = c(xlim, 5), yaxp = c(ylim, 5))
#    axis(1)
#    axis(2)
#}

#paper.plot <- function(m, ...) {
#  plot(m, initplot=paper.initplot, transform=dlp.transform, pch=20, ...)
#}

###

kld1 <- function(m, tr) {
  nk <- length(m$active)
  km <- matrix(0, nk, nk)
  for (i in 1:nk) {
    k <- m$active[i]
    tk <- tr$tronly[m$Z==k,1]
    pk <- sum(tk==0)/length(tk)
    for (j in 1:nk) {
      l <- m$active[j]
      tl <- tr$tronly[m$Z==l,1]
      pl <- sum(tl==0)/length(tl)
      km[i,j] <- pk*log(pk/pl) + (1-pk)*log((1-pk)/(1-pl))
    }
  }
  return(km)
}
