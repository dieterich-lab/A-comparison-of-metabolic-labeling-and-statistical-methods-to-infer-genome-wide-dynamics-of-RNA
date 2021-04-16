#!/biosw/R/3.5.1/bin/Rscript

# Usage: ./get_estimates.R 

# Compute decay rates and confidence intervals for GRAND-SLAM for the different protocols.
# MAP estimated using matched subset of time points (pulseR), approximate CIs obtained 
# for each subset using the calculated MAP

library(dplyr)
library(tibble)
library(purrr)

loc <- here::here()

# GRAND-SLAM results
gsLoc <- file.path(loc, "grand-slam", fsep=.Platform$file.sep)
gsFit <- "all"

# pulseR results 
pulseDir <- file.path(loc, "pulseRTc", "results", fsep=.Platform$file.sep)
pulseFit <- "pulsefit"
pulseCis <- "pulsecis"

# Output
tabDir <- file.path(loc, "paper", "tables", fsep=.Platform$file.sep)
gsDir <- file.path(gsLoc, "tables", fsep=.Platform$file.sep)

# Fitting sets - 3 points or all
prefix <- "../pulser"
source(file.path(loc, "pulseRTc", "pulser", "utils.R", fsep=.Platform$file.sep)) 
timeSets <- timeSets3
timeSets[length(timeSets)+1] <- allSets

dataSets <- c("ercc", "std", "slam", "tls", "tuc")

## functions

# log likelihood
loglik <- function(d, pars) sum((pars$a-1)*log(1-exp(-pars$t*d))-pars$t*d*pars$b)

# approximate confidence intervals
getCI <- function(ll, interval, pars, optimum, confidence=0.95) {
  threshold <- stats::qchisq(confidence, 1)/2
  objective <- function(x) ll(optimum, pars) - ll(x, pars) - threshold
  optimalObjective <- objective(optimum)
  ci <- c(NA, NA)
  if (optimalObjective * objective(interval[1]) < 0)
    ci[1] <- stats::uniroot(objective, c(interval[1], optimum))$root
  if (optimalObjective * objective(interval[2]) < 0)
    ci[2] <- stats::uniroot(objective, c(optimum, interval[2]))$root
  ci
}

# maximum a posteriori estimate for decay rate
mle <- function(v, use=1:8, interval=c(1e-12,2), confidence=0.95) {
    # fixed output format GRAND-SLAM
    a <- v[seq(1, length(v), by=2)]
    b <- v[seq(2, length(v), by=2)]
    # all time points for one given param
    t <- c(1, 2, 4, 8, 1, 2, 4, 8)
    pars <- list()
    pars$a <- a[grep(use, names(a))]
    pars$b <- b[grep(use, names(b))]
    pars$t <- t[grep(use, names(a))]
    if(any(is.nan(c(pars$a,pars$b)))) return(list(map=NA, ci=c(NA, NA)))
    optimum <- optimize(loglik, interval, pars, maximum=T)$maximum
    ci <- getCI(loglik, interval, pars, optimum, confidence)
    list(map=optimum, ci=ci)
}


## data

gsFiles <- file.path(gsLoc , dataSets[grep("slam|tuc|tls", dataSets)], gsFit, paste(gsFit, "tsv", sep="."), fsep=.Platform$file.sep)
print(paste("Using: ", gsFiles, " ...", sep=""))

gsList <- map(gsFiles, function(.id) {
    # only read alpha/beta
    gs <- read.table(.id, sep='\t', header=T, row.names = 1)
    gs[,grepl('alpha|beta', colnames(gs))]
})

# select for each data the matching time sets
MAPr <- map(gsList, 
            function(.id) {
                ld <- list()
                data <- unlist(strsplit(colnames(.id)[1], "_"))[2]
                for (idx in seq_along(timeSets)) {
                    d <- data.frame(gene=rownames(.id), check.names=FALSE)
                    tp <-  as.integer(timeSets[[idx]])
                    tp <- tp[!tp==0] # not fitted by GS, used to model pe
                    label <- paste(data, paste(paste(tp, collapse=","), "h", sep=""))
                    use <- paste(paste(tp, 'h', sep = ""), collapse = "|")
                    out <- unlist(apply(as.matrix(.id), 1, mle, use=use), 
                                  recursive = F, use.names = F) 
                    cis <- unlist(out[seq(2, length(out), by=2)])
                    d[[paste(paste('GS', label), "d", sep = " ")]] <- as.numeric(out[seq(1, length(out), by=2)])
                    d[[paste(paste('GS', label), "d.min", sep = " ")]] <- cis[seq(1, length(cis), by=2)]
                    d[[paste(paste('GS', label), "d.max", sep = " ")]] <- cis[seq(2, length(cis), by=2)]
                    rownames(d) <- d$gene
                    d$gene <- NULL
                    saveRDS(d,
                        sprintf(file.path(gsDir, paste(tolower(data), "-%s.rds", sep="")),
                        paste(tp, collapse ="-")))
                    ld[[label]] <- d
                }
                ld
})
gsList <- unlist(MAPr, recursive = FALSE)

# now we also combine these results with the pulseR estimates and generate a table
# with all common genes across all experiments for each time set

# first get pulseR degradation rates
pulseFiles <- lapply(timeSets, function(x) paste(paste(pulseFit, paste(x, collapse = "-"), sep = "-"), "rds", sep = "."))
pulseFiles <- file.path(pulseDir , rep(dataSets, each = length(pulseFiles)), pulseFiles, fsep=.Platform$file.sep)

print(paste("Using: ", pulseFiles, " ...", sep=""))

pulseList <- map(pulseFiles, function(.id) {
  data <- toupper(rev(unlist(strsplit(.id, "/")))[2])
  # reformat
  label <- paste(data, paste(paste(paste(unlist(strsplit(gsub(".rds", "", basename(.id)), '-'))[-1], collapse=","), collapse = " "), "h", sep=""))
  pulse <- readRDS(.id)
  pulsedf <- as.data.frame(pulse$fit$d)
  cis <- readRDS(gsub(pulseFit, pulseCis, .id))
  cis <- as.data.frame(cis)
  df <- cbind(pulsedf, cis)
  colnames(df) <- c(paste(paste("pulseR", label), "d", sep = " "), 
                    paste(paste("pulseR", label), "d.min", sep = " "), 
                    paste(paste("pulseR", label), "d.max", sep = " "))
  rownames(df) <- rownames(pulse$pd$counts)
  df
})

pulseNames <- map(pulseFiles, function(.id) {
  data <- toupper(rev(unlist(strsplit(.id, "/")))[2])
  # reformat
  paste(data, paste(paste(paste(unlist(strsplit(gsub(".rds", "", basename(.id)), '-'))[-1], collapse=","), collapse = " "), "h", sep=""))
})
names(pulseList) <- pulseNames

# match everything and write to disk
res <- lapply(timeSets, 
              function(ts) {
              # pulseR
              tp <- paste(paste(ts, collapse=","), "h", sep="")
              p <- pulseList[grep(paste("^.{1,5}", tp, "$", sep = ""), names(pulseList))]
              # GRAND-SLAM
              tp <- ts[!ts==0]
              tp <- paste(paste(tp, collapse=","), "h", sep="")
              g <- gsList[grep(paste("^.{1,5}", tp, "$", sep = ""), names(gsList))]
              matched <- c(p, g)
              # first remove where d is NA (GS)    
              matched <- map(matched, function(df) { 
                      col <- colnames(df)[grep('min|max', colnames(df), invert=T)]
                      df[!is.na(df[[col]]),]
              })
              genes <- map(matched, function(df) { rownames(df) })
              idx <- Reduce(intersect, genes)
              matched <- map(matched, function(df) { df[match(idx, rownames(df)),] })
              names(matched) <- NULL
              matched_df <- do.call(cbind, matched)
              saveRDS(matched_df,
                sprintf(file.path(tabDir, "tbl-%s.rds"), paste(ts, collapse ="-")))
              matched_df
})

resNames <- lapply(timeSets, function(ts) { paste(ts, collapse="-") })
names(res) <- resNames

             
# rewrite as xlsx
library(openxlsx)

workBook <- createWorkbook()
w <- lapply(seq_along(res),
            function(.id) {
            n <- names(res)[.id]
            r <- res[[.id]] 
            addWorksheet(workBook, sheetName=n)
            writeDataTable(workBook, sheet=.id, x=r, rowNames = TRUE)
})
saveWorkbook(workBook, file.path(tabDir, "tbl.xlsx"), overwrite=FALSE)
