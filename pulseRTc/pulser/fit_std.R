#! /usr/bin/env Rscript

# Usage: ./fit_std.R <0/1> <0/1>
# <0/1>: 0 = fit (default), 1 = use ERCC spike-ins
# <0/1>: 0 = all time points (default), 1 = subset of time points

library(pulseR)

library(purrr)
library(tibble)
library(dplyr)


## Adapted from https://github.com/dieterich-lab/DesignMetabolicRNAlabeling
## Fit the kinetic model to the read counts including total fraction at 0h, 
## with or without ERCC spike-ins
## pulseR_1.0.3

loc <- here::here("pulseRTc", "pulser")

source(file.path(loc, "utils.R", fsep=.Platform$file.sep)) 
source(file.path(loc, "models.R", fsep=.Platform$file.sep))

## local options
expressionThreshold <- 50

## input/output: all paths are relative to this directory

whichData <- "std"

args <- commandArgs(trailingOnly=TRUE)
use_spikes <- FALSE
use_fit <- FALSE
usedSets <- allSets
modelStr <- "pulse"
if (length(args)==1 & as.integer(args[1])==1) { 
    use_spikes <- TRUE
    whichData <- "ercc"
} else if (length(args)>=1) {
    if (as.integer(args[1])==1) {
        use_spikes <- TRUE
        whichData <- "ercc"
    }
    if (as.integer(args[2])==1) {
        # we can use sets with 2 points...
        # usedSets <- timeSets2
        usedSets <- timeSets3
        use_fit <- TRUE
    }
}


prefix <- '../results'
rdsDir <- file.path(loc, prefix, whichData, fsep=.Platform$file.sep)
print(paste("Fitting ", rdsDir, " ...", sep=""))

## data

# make sure "std" count data is linked under ercc if using spike-ins... 
data <- getPulseData(rdsDir, refLevel="total")
counts <- data$counts
conditions <- data$conditions


## Filter out lowly expressed genes (incl. ERCC spike-ins)
## based on total counts (total at 0h)

highExpr <- whichHigh(1 + counts[, conditions$fraction == "total"],
  expressionThreshold)
counts <- counts[highExpr,]


## fitting options: tolerance and upper/lower bounds for parameters

tolerance <- list(params      = 0.01,
                  normFactors = 0.01,
                  logLik      = 0.01)
                
boundaries <- list(mu          = log(c(1e-2, 1e6)),
                   d           = c(1e-3, 2),
                   size        = c(1, 1e3),
                   normFactors = c(1e-6, 20))
                       

if (use_spikes) {
    ## Prepare spike-ins list...
    
    erccList <- rownames(counts[grep('^ERCC-', rownames(counts)),])
    print(paste("Using ", length(erccList) , " ERCC spike-ins ...", sep=""))
    spikeIns <- list() 
    spikeIns$refGroup <- "total"
    spikeList <- list()
    spikeList$total <- list(erccList)
    spikeList$flow_through <- list(erccList)
    spikeList$pull_down <- list(erccList)
    spikeIns$spikeLists <- spikeList
    
    
    ## adjust
    
    tolerance$normFactors <- NULL
    boundaries$normFactors <- NULL
    
    
    ## define function to set initial values
    # with spike-ins, normFactors are not created
    init <- function(pd, opts) {
        fit <- initParameters(list(),
                        geneParams = c("mu", "d"),
                        pulseData = pd,
                        opts)
        fit$mu <- log(1e-1 + pd$counts[,1])
        fit$size <- 1e2
        ## use prior fit
        if (use_fit) {
            fit_str <- sprintf(paste(modelStr, "fit-%s.rds", sep=""), do.call("paste", c(allSets, collapse="-")))
            fit_str <- file.path(rdsDir, fit_str, fsep=.Platform$file.sep)
            print(paste("Using ", fit_str, " as initial conditions...", sep=""))
            fit <- readRDS(fit_str)$fit
        }
        fit
    }

    
} else {
    ## or remove spike-ins
    
    print("Fitting normFactors ... ")
    counts <- counts[!grepl("^ERCC-", rownames(counts)),]
    spikeIns <- NULL

    
    ## define function to set initial values
    init <- function(pd, opts) {
        fit <- initParameters(list(),
                        geneParams = c("mu", "d"),
                        pulseData = pd,
                        opts)
        fit$mu <- log(1e-1 + pd$counts[,1])
        fit$size <- 1e2
        # initialize normFactors
        fit$normFactors <-  lapply(fit$normFactors, function(x) {x[1] <- 1; x})
        ## use prior fit
        if (use_fit) {
            fit_str <- sprintf(paste(modelStr, "fit-%s.rds", sep=""), do.call("paste", c(allSets, collapse="-")))
            fit_str <- file.path(rdsDir, fit_str, fsep=.Platform$file.sep)
            print(paste("Using ", fit_str, " as initial conditions...", sep=""))
            fit <- readRDS(fit_str)$fit
        }
        fit
    }
    
}


## model formulas

fractions <- c("total", "flow_through", "pull_down")
model <- makeFormStdPulse(fractions)


## run pulseR

getFit <- function(counts, conditions, model, init, norms, boundaries, tolerance, spikes) {
  pd <- makePD(counts, conditions, model, norms=norms, spikes=spikes)
  opts <- setOpts(boundaries, tolerance, normFactors=boundaries$normFactors)
  initf <- match.fun(init)
  initPars <- initf(pd, opts) # use pulseData here
  fit <- fitModel(pd, initPars, opts)
  list(fit=fit, pd=pd, opts=opts)
}


fitTimePoints <- function(counts, conditions, model, init, norms, boundaries, tolerance, spikes, tp) {
  tpconditions <- conditions[conditions$time %in% tp,]
  tpcounts <- counts[, tpconditions$sample]
  res <- getFit(tpcounts, tpconditions, model, init, norms, boundaries, tolerance, spikes)
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
  norms=NULL, # for rt-con data 
  boundaries=boundaries, 
  tolerance=tolerance,
  spikes=spikeIns)

  
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
  
