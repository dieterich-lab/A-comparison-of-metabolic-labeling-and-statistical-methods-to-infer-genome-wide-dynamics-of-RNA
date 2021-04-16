#!/biosw/R/3.5.1/bin/Rscript

library(purrr)
library(tibble)
library(dplyr)


##' Prepare count tables and metadata for analysis
##' Note: this function is not all purpose, we know the exact format of the data!
##' No filtering based on expression, etc.
getData <- function(data, dataDir, dataFile, outDir){
    cts <- read.table(file.path(dataDir, dataFile, fsep=.Platform$file.sep), row.names=1, sep = "\t", header=TRUE, check.names=FALSE)
    cts <- cts %>% dplyr::select(starts_with("/prj"))
    depth <- length(strsplit(colnames(cts)[1], "/", fixed=T)[[1]])
    conditions <- data.frame(vapply(strsplit(colnames(cts), "/", fixed=T), "[", "", depth), stringsAsFactors=FALSE)
    colnames(conditions) <- "sample"
    if (data == 'std') {
        conditions$fraction <- vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 2)
        conditions$fraction[grepl("flow", conditions$fraction)] <- "flow_through"
        conditions$fraction[grepl("enriched", conditions$fraction)] <- "pull_down"
        conditions$time <- as.numeric(gsub("h", "", gsub("^.*_", "", vapply(strsplit(conditions$sample, ".", fixed=T), "[", "", 1))))
        conditions$rep <- gsub("^[1-9][0-9][0-9][0-9][0-9][0-9]", "", vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 1))
        conditions$sample <- gsub("[A-Z]$", "", vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 1))
        # redefine columns
        colnames(cts) <- conditions$sample
        # remove spike-in mix, but keep ERCC
        # we know which ones
        spikes <- c("CAT", "Luc", "Rluc_with-4-point-mutations")
        cts <- cts[!rownames(cts) %in% spikes, ]
        # if we use fractional counts (multimappers)
        cts <- round(cts)
    } else { # slam, tls, or tuc
        conditions$fraction <- vapply(strsplit(vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 3), ".", fixed=T), "[", "", 2)
        conditions$time <- as.numeric(gsub("h", "", vapply(strsplit(vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 3), ".", fixed=T), "[", "", 1)))
        conditions$rep <- gsub("^[1-9][0-9][0-9][0-9][0-9][0-9]", "", vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 1))
        conditions$sample <- gsub("[A-Z]$", "", vapply(strsplit(conditions$sample, "_", fixed=T), "[", "", 1))
        conditions$sample <- paste(conditions$sample, "L", sep="")
        conditions$sample[grepl("unlabelled", conditions$fraction)] <- gsub("L", "U", conditions$sample[grepl("unlabelled", conditions$fraction)])
        # redefine columns
        colnames(cts) <- conditions$sample
        # remove spike-in mix and ERCC (unused)
        # we know which ones
        spikes <- c("CAT", "Luc", "Rluc_with-4-point-mutations")
        cts <- cts[!rownames(cts) %in% spikes, ]
        cts <- cts[!grepl("ERCC-", rownames(cts)),]
        # if we use fractional counts (multimappers)
        cts <- round(cts)
    }
    saveRDS(cts, file.path(outDir, "counts.rds", fsep=.Platform$file.sep))
    saveRDS(conditions, file.path(outDir, "samples.rds", fsep=.Platform$file.sep))
}


## Call

loc <- here::here("pulseRTc", "pulser")
dataDir <- '../workflow/tables'
dataDir <- file.path(loc, dataDir, fsep=.Platform$file.sep)
print(paste("Processing ", dataDir, " ...", sep=""))

outDir <- '../results'
outDir <- file.path(loc, outDir, fsep=.Platform$file.sep)

tsv <- "-0h-1h-2h-4h-8h.counts.tsv"

lapply(c("std", "slam", "tuc", "tls"), 
        function(data) getData(data, dataDir, paste(data, tsv, sep=""), file.path(outDir, data, "data", fsep=.Platform$file.sep)))


