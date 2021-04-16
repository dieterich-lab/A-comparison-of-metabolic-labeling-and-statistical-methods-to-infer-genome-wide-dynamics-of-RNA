## Adapted from https://github.com/dieterich-lab/DesignMetabolicRNAlabeling


##' Construct all rds result files
getPath <- function(rdsDir) {
  x <- c("pd", "fit", "opts", "cis")
  path <- lapply(setNames(x, x),
                    function(y)
                      file.path(rdsDir, paste0(y, ".rds")))
  path
}


##' Read counts (counts.rds) and metadata (samples.rds)
##' Set reference factor level
getPulseData <- function(rdsDir, refLevel=NULL) {
  cts <- readRDS(file.path(rdsDir, "data", "counts.rds"))
  conditions <- readRDS(file.path(rdsDir, "data", "samples.rds"))
  # ordering columns to match metadata
  cts <- cts[conditions$sample]
  conditions <- conditions[c("fraction", "time", "sample", "rep")]
  conditions$fraction <- factor(conditions$fraction)
  if (!is.null(refLevel)) { 
    conditions$fraction <-
        relevel(conditions$fraction, refLevel)
  }
  list(conditions = conditions,
       counts = cts)
}


##' Filter out lowly expressed genes using geometric mean
whichHigh <- function(x, level) {
  apply(x, 1, function(y) exp(mean(log(y)))) > level
}


##' Set all options
setOpts <- function(bounds, tolerance, cores = 24, replicates = 5, normFactors = NULL) {
  opts <- setFittingOptions(verbose = "verbose")
  opts$cores <- cores
  opts$replicates <- replicates
  ## if rt-conversion data, we do not fit normalisation coefficients, because they
  ## are derived from the DESeq-like normalisation
  if (is.null(normFactors)) { opts$fixedNorms <- TRUE }
  opts <- setBoundaries(bounds, normFactors = normFactors, options = opts)
  opts <- setTolerance(
    params = tolerance$params,
    normFactors = tolerance$normFactors,
    logLik = tolerance$logLik,
    options = opts
  )
  if (is.null(tolerance$normFactors)) {
    opts$tolerance$normFactors <- NULL
  }
  opts
}


##' Create PulseData object
makePD <- function(counts, conditions, formulas, norms=NULL, spikes=NULL) {
  if (!is.null(spikes)) {
    pd <- PulseData(
        counts[, conditions$sample],
        conditions[c("fraction", "time")], #  data.frame; the first column corresponds to the conditions given in formulas
        formulas$formulas,
        formulas$formulaIndexes,
        spikeins=spikes # passed, and removed from count table
    )
  } else {
    pd <- PulseData(
        counts[, conditions$sample],
        conditions[c("fraction", "time")], #  data.frame; the first column corresponds to the conditions given in formulas
        formulas$formulas,
        formulas$formulaIndexes,
        groups=~fraction+time
    )
  }
  ## if rt-conversion data, add normalisation factors
  if (!is.null(norms)) {
    print("Adding normalisation factors for rt-conversion data...")
    pd$depthNormalisation <- norms[conditions$sample]
  }
  pd
}


##' Fitting sets
# for RTC data
timeSets3 <- list(
  c(0, 1, 2),
  c(0, 2, 4),
  c(0, 4, 8)
)
# for standard
timeSets2 <- list(
  c(0, 1),
  c(0, 2),
  c(0, 4),
  c(0, 8)
)
allSets <- list(
  c(0, 1, 2, 4, 8)
)


##' a helper to read fit results and to make a result table.
##' returns a list with
##'     result$pd  - PulseData object
##'     result$fit - fitting results
##'     result$tab - fitting results with confidence intervals as a data.frame
readFits <- function(d) {
  res <- readRDS(file.path(d, "fit.rds"))
  pd <- readRDS(file.path(d, "pd.rds"))
  fits <- res[c("mu", "d")]
  if (file.exists(file.path(d, "cis.rds"))) {
    cis <- readRDS(file.path(d, "cis.rds"))
    cis <- lapply(cis,
                  function(x) {
                    setNames(as.data.frame(x), c("min", "max"))
                  })
    cis <- do.call(cbind, cis)
    fits <- cbind(as.data.frame(fits), cis)
  }
  list(fit = res,
       pd = pd,
       tab = as.data.frame(fits))
}


##' Create a data frame with fitted params and their CI
sumup <- function(fit, cis, pd = NULL) {
  fit$normFactors <- NULL
  fit$size <- NULL
  fit <- fit[vapply(fit, length, integer(1)) > 1]
  fit <- as.data.frame(fit)
  fit$id <- rownames(pd$counts)
  cisNames <-
    lapply(names(cis), function(x)
      paste(x, c("min", "max"), sep = "."))
  cis <- setNames(as.data.frame(do.call(cbind, cis)), unlist(cisNames))
  cbind(fit, cis)
}

