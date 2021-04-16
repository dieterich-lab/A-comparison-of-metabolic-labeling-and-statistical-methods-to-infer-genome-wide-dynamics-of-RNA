#!/biosw/R/3.5.1/bin/Rscript

# Usage: ./fit_rtc.R [slam/tuc/tls] <0/1>
# <0/1>: 0 = all time points (default), 1 = subset of time points

library(pulseR)

library(purrr)
library(tibble)
library(dplyr)

## Adapted from https://github.com/dieterich-lab/DesignMetabolicRNAlabeling
## Fit the kinetic model to the read counts for reverse transcription nucleotide 
## conversion data: SLAM, TUC, TLS.
## pulseR_1.0.3

loc <- here::here("pulseRTc", "pulser")

source(file.path(loc, "utils.R", fsep=.Platform$file.sep)) 
source(file.path(loc, "models.R", fsep=.Platform$file.sep))

## local options
expressionThreshold <- 50


## input/output: all paths are relative to this directory

args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1) { stop("./fit_rtc.R [slam/tuc/tls] <0/1>\n", call.=FALSE) }
whichData <- tolower(args[1])

use_fit <- FALSE
usedSets <- allSets
modelStr <- "pulse"
if (length(args)>1 & as.integer(args[2])==1) { 
    usedSets <- timeSets3
    use_fit <- TRUE
}


prefix <- '../results'
rdsDir <- file.path(loc, prefix, whichData, fsep=.Platform$file.sep)
print(paste("Fitting ", rdsDir, " ...", sep=""))

## data

data <- getPulseData(rdsDir)
counts <- data$counts
conditions <- data$conditions


## Filter out lowly expressed genes, and determine normalisation factors

conditions$map <- gsub('L|U', '', conditions$sample)
totals <- map(unique(conditions$map), function(.id) {
  i <- conditions$sample[conditions$map == .id]
  apply(counts[,i], 1, sum)
})
totals <- do.call(cbind, totals)

highExpr <- whichHigh(1 + totals, expressionThreshold)
counts <- counts[highExpr,]
totals <- totals[highExpr,]

norms <- pulseR:::findDeseqFactorsSingle(totals)
names(norms) <- unique(conditions$map)
norms <- norms[conditions$map]
names(norms) <- conditions$sample


## fitting options: tolerance and upper/lower bounds for parameters

tolerance <- list(params = 0.01,
                  logLik = 0.01)

# assume 20:1 labelling efficiency, neglecting substitution, ratio is approx. the same
boundaries <- list(mu1  = c(log(1e-2*20), log(1e6*20)), # substituted variable
                   mu2  = c(log(1e-2), log(1e6)),
                   mu3  = c(log(1e-2), log(1e6)),
                   d    = c(1e-3, 2),
                   size = c(1, 1e3))
                   
## set initial values here 

init <- function(counts) {
  fit <- list(
    mu1  = log(1e-1 + counts[,1]),
    mu2  = log(1e-1 + counts[,1]),
    mu3  = log(1e-1 + counts[,1]),
    d    = rep(5e-1, length(counts[,1])),
    size = 1e2)
  ## use prior fit
  if (use_fit) {
        fit_str <- sprintf(paste(modelStr, "fit-%s.rds", sep=""), do.call("paste", c(allSets, collapse="-")))
        fit_str <- file.path(rdsDir, fit_str, fsep=.Platform$file.sep)
        print(paste("Using ", fit_str, " as initial conditions...", sep=""))
        fit <- readRDS(fit_str)$fit
    }
  fit
}


## model formulas

model <- makeFormRtcPulse()

## run pulseR

getFit <- function(counts, conditions, model, init, norms, boundaries, tolerance) {
  pd <- makePD(counts, conditions, model, norms=norms)
  opts <- setOpts(boundaries, tolerance)
  initf <- match.fun(init)
  initPars <- initf(pd$counts) # use pulseData here
  fit <- fitModel(pd, initPars, opts)
  list(fit=fit, pd=pd, opts=opts)
}


fitTimePoints <- function(counts, conditions, model, init, norms, boundaries, tolerance, tp) {
  tpconditions <- conditions[conditions$time %in% tp,]
  tpcounts <- counts[, tpconditions$sample]
  res <- getFit(tpcounts, tpconditions, model, init, norms, boundaries, tolerance)
  saveRDS(
    res,
    sprintf(file.path(rdsDir, paste(modelStr, "fit-%s.rds", sep="")), 
        paste(tp, collapse ="-")))
  res
}

# Adjust boundaries and init depending on model

res <- map(
  usedSets,
  fitTimePoints,
  counts=counts, 
  conditions=conditions, 
  model=model, 
  init=init,
  norms=norms, 
  boundaries=boundaries, 
  tolerance=tolerance)

## compute confidence intervals for d

cisr <- map(
  res,
  function(.x) {
    cis <- ciGene("d",
      par = .x$fit,
      geneIndexes = seq_along(.x$fit$d),
      pd = .x$pd,
      options = .x$opts)
    tp <- unique(.x$pd$conditions$time)
    saveRDS(cis,
      sprintf(file.path(rdsDir, paste(modelStr, "cis-%s.rds", sep="")),
        paste(tp, collapse ="-")))
    cis
  })
  
