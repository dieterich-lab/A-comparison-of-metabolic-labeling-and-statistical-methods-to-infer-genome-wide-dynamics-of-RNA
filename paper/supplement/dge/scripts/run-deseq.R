#!/biosw/R/3.6.2/bin/Rscript

## Usage: ./run-deseq.R [data]
## data: slam, tls, or tuc

library(DESeq2)
library(IHW)

library("Glimma")
library("genefilter")
library("org.Hs.eg.db")

library(dplyr)
library(purrr)
library(tibble)

library(openxlsx)

# ---------------------------------------------------------

## Default LFC and FDR threshold

lfcThreshold.set <- log2(1.2)
altHypothesis.set <- "greaterAbs"
alpha.set <- 0.05

# ---------------------------------------------------------

## functions

# filter out lowly expressed genes using geometric mean
whichHigh <- function(x, level) {
  apply(x, 1, function(y) exp(mean(log(y)))) > level
}

# DESeq results
get_results <- function (contrast, dds) { 

    num <- as.character(contrast[1])
    denom <- as.character(contrast[2]) # reference 

    print(paste("Contrast: time_", num, "_vs_", denom, sep=""))

    res <- results(dds,
                   contrast=c("time", num, denom),
                   lfcThreshold=lfcThreshold.set, 
                   altHypothesis=altHypothesis.set,
                   alpha=alpha.set,
                   filterFun=ihw)
    res$padj[is.na(res$padj)] <- 1
                   
    res$symbol <- mapIds(org.Hs.eg.db,
                         keys=rownames(res),
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")

    res.shrunken <- lfcShrink(dds, 
                              coef=resultsNames(dds)[grepl(paste('time_', num, sep=''), resultsNames(dds))], 
                              res=res,
                              type="apeglm", 
                              lfcThreshold=lfcThreshold.set)
    res.shrunken$pvalue <- res$pvalue
    res.shrunken$padj <- res$padj
    
    res.shrunken$symbol <- mapIds(org.Hs.eg.db,
                                  keys=rownames(res.shrunken),
                                  column="SYMBOL",
                                  keytype="ENSEMBL",
                                  multiVals="first")
                         
    is.de <- as.numeric(res.shrunken$padj < alpha.set & abs(res.shrunken$log2FoldChange) > lfcThreshold.set)
    anno <- data.frame(GeneID=rownames(res.shrunken), symbol=res.shrunken$symbol)
    glMDPlot(res.shrunken, 
             counts=counts(dds ,normalized=TRUE),
             anno,
             dds$time, 
             samples=colnames(dds), 
             status=is.de, 
             transform = FALSE,
             xlab = "logMeanExpr",
             ylab = "log2FoldChange",
             side.ylab = "NormalizedCount",
             path=dirloc.out, 
             folder=paste("glimma-plots", num, "_vs_", denom, sep=""), 
             launch=FALSE)

    res.tib <- res %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>% 
        as_tibble()
    
    res.shrunken.tib <- res.shrunken %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>% 
        as_tibble()
        
    # write to disk
    wb <- createWorkbook()

    addWorksheet(wb, sheetName=paste0(num, "_vs_", denom, sep=""))
    writeDataTable(wb, sheet=1, x=res.tib)

    addWorksheet(wb, sheetName=paste0(num, "_vs_", denom, "_shrunken", sep=""))
    writeDataTable(wb, sheet=2, x=res.shrunken.tib)
    
    filen <- paste0("time_", num, "_vs_", denom, ".xlsx", sep="")
    filen <- file.path(dirloc.out, filen, fsep=.Platform$file.sep)
    saveWorkbook(wb, filen, overwrite=TRUE)
    
}


# ---------------------------------------------------------

## Call

args <- commandArgs(trailingOnly=TRUE)
data <- args[1]

loc <- here::here("dge")
dataDir <- file.path(loc, 'tables')
tsv <- "-0h-1h-2h-4h-8h.counts.tsv"
dataFile <- paste(data, tsv, sep="")
print(paste("Processing ", dataFile, " ...", sep=""))

dirloc.out <- file.path(loc, 'results', data)
dir.create(dirloc.out, showWarnings=FALSE)
print(paste("Writing to ", dirloc.out, " ...", sep=""))


cts <- read.table(file.path(dataDir, dataFile, fsep=.Platform$file.sep), row.names=1, sep = "\t", header=TRUE, check.names=FALSE)
cts <- cts %>% dplyr::select(starts_with("/prj"))
depth <- length(strsplit(colnames(cts)[1], "/", fixed=T)[[1]])
coldata <- data.frame(vapply(strsplit(colnames(cts), "/", fixed=T), "[", "", depth), stringsAsFactors=FALSE)
colnames(coldata) <- "sample"

# here we're not using time as a continuous variable e.g. to identify genes that change as a function of time
# but rather want to compare the time points in terms of expression, so we use them as factors/contrasts
# coldata$time <- as.numeric(gsub("h", "", vapply(strsplit(vapply(strsplit(coldata$sample, "_", fixed=T), "[", "", 3), ".", fixed=T), "[", "", 1)))

if (data == 'std') {
    coldata$fraction <- vapply(strsplit(coldata$sample, "_", fixed=T), "[", "", 2)
    coldata$fraction[grepl("flow", coldata$fraction)] <- "flow_through"
    coldata$fraction[grepl("enriched", coldata$fraction)] <- "pull_down"
    coldata$time <- gsub("h", "", gsub("^.*_", "", vapply(strsplit(coldata$sample, ".", fixed=T), "[", "", 1)))
    coldata$time <- as.factor(coldata$time)
    coldata$rep <- gsub("^[1-9][0-9][0-9][0-9][0-9][0-9]", "", vapply(strsplit(coldata$sample, "_", fixed=T), "[", "", 1))
    coldata$rep <- as.factor(coldata$rep)
    coldata$sample <- gsub("[A-Z]$", "", vapply(strsplit(coldata$sample, "_", fixed=T), "[", "", 1))
    rownames(coldata) <- coldata$sample
    colnames(cts) <- coldata$sample

    # remove 1h time point, missing from IP A
    coldata <- coldata[!coldata$time=='1',]
    coldata$time <- droplevels(coldata$time)
    cts <- cts[,colnames(cts) %in% rownames(coldata)]

    # use "total samples"
    total <- cts[,colnames(cts) %in% coldata[coldata$time=='0',]$sample]
    colnames(total) <- c('0A', '0B')

    reduced <- coldata[!coldata$time=='0',]
    reduced$time <- droplevels(reduced$time)
    reduced$map <- paste(reduced$time, reduced$rep, sep='')
    cts <- map(unique(reduced$map), function(.id) {
    i <- reduced$sample[reduced$map == .id]
    apply(cts[,i], 1, sum)
    })

    cts <- do.call(cbind, cts)
    cts <- cbind(cts, total)
    cts <- cts[,c(7,1,2,3,8,4,5,6)]

    # now adjust coldata
    coldata$sample <- paste(coldata$time, coldata$rep, sep='')
    coldata <- coldata[!duplicated(coldata$sample),]
    coldata$time <- droplevels(coldata$time)
    coldata$rep <- droplevels(coldata$rep)
    coldata$fraction <- NULL
    rownames(coldata) <- coldata$sample

    colnames(cts) <- coldata$sample
} else {
    coldata$time <- gsub("h", "", vapply(strsplit(vapply(strsplit(coldata$sample, "_", fixed=T), "[", "", 3), ".", fixed=T), "[", "", 1))
    coldata$time <- as.factor(coldata$time)
    coldata$rep <- gsub("^[1-9][0-9][0-9][0-9][0-9][0-9]", "", vapply(strsplit(coldata$sample, "_", fixed=T), "[", "", 1))
    coldata$rep <- as.factor(coldata$rep)
    coldata$sample <- gsub("[A-Z]$", "", vapply(strsplit(coldata$sample, "_", fixed=T), "[", "", 1))
    rownames(coldata) <- coldata$sample
    colnames(cts) <- coldata$sample
}

# remove spike-in mix and ERCC (unused)
# we know which ones
spikes <- c("CAT", "Luc", "Rluc_with-4-point-mutations")
cts <- cts[!rownames(cts) %in% spikes, ]
cts <- cts[!grepl("ERCC-", rownames(cts)),]
# if we use fractional counts (multimappers)
cts <- round(cts)

# expression threshold - same as for pulseR
# note however that the `results` function performs independent filtering by default using the
# mean of normalized counts as a filter statistic
expressionThreshold <- 50
highExpr <- whichHigh(1 + cts, expressionThreshold)
cts <- cts[highExpr,]

stopifnot(all(rownames(coldata) == colnames(cts)))

# contrasts to test
contrasts <- data.frame(c('1', '2', '4', '8'), c('0', '0', '0', '0'))
colnames(contrasts) <- c('cond', 'ref')

if (data == 'std') {
    contrasts <- contrasts[!contrasts$cond=='1',]
}

# construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=~time)

# fit all at once
dds$time <- relevel(dds$time, ref='0')
dds <- DESeq(dds)
                                          
apply(contrasts, 1, get_results, dds=dds)
