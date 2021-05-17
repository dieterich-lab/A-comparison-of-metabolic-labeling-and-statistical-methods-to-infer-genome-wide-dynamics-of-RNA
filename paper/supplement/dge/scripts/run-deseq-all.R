#! /usr/bin/env Rscript

## Usage: ./run-deseq-all.R

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

getData <- function(data, dataDir, dataFile) {
    print(paste("Processing ", dataFile, " ...", sep=""))
    cts <- read.table(file.path(dataDir, dataFile, fsep=.Platform$file.sep), row.names=1, sep = "\t", header=TRUE, check.names=FALSE)
    cts <- cts %>% dplyr::select(starts_with("/prj"))
    depth <- length(strsplit(colnames(cts)[1], "/", fixed=T)[[1]])
    coldata <- data.frame(vapply(strsplit(colnames(cts), "/", fixed=T), "[", "", depth), stringsAsFactors=FALSE)
    colnames(coldata) <- "sample"
    coldata$time <- gsub("h", "", vapply(strsplit(vapply(strsplit(coldata$sample, "_", fixed=T), "[", "", 3), ".", fixed=T), "[", "", 1))
    coldata$rep <- gsub("^[1-9][0-9][0-9][0-9][0-9][0-9]", "", vapply(strsplit(coldata$sample, "_", fixed=T), "[", "", 1))
    coldata$sample <- gsub("[A-Z]$", "", vapply(strsplit(coldata$sample, "_", fixed=T), "[", "", 1))
    coldata$data <- data
    rownames(coldata) <- coldata$sample
    colnames(cts) <- coldata$sample
    spikes <- c("CAT", "Luc", "Rluc_with-4-point-mutations")
    cts <- cts[!rownames(cts) %in% spikes, ]
    cts <- cts[!grepl("ERCC-", rownames(cts)),]
    # if we use fractional counts (multimappers)
    cts <- round(cts)
    list("coldata"=coldata, "cts"=cts)
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

loc <- here::here("dge")
dirloc.out <- file.path(loc, 'results', 'all')
dir.create(dirloc.out, showWarnings=FALSE)

dataDir <- file.path(loc, 'tables')
tsv <- "-0h-1h-2h-4h-8h.counts.tsv"

l <- lapply(c("slam", "tuc", "tls"), 
        function(data) getData(data, dataDir, paste(data, tsv, sep="")))

coldata <- do.call(rbind, lapply(l, `[[`, 1))
coldata$time <- as.factor(coldata$time)
coldata$rep <- as.factor(coldata$rep)
coldata$data <- as.factor(coldata$data)

cts <- do.call(cbind, lapply(l, `[[`, 2))

stopifnot(all(rownames(coldata) == colnames(cts)))

print(paste("Writing to ", dirloc.out, " ...", sep=""))

expressionThreshold <- 50
highExpr <- whichHigh(1 + cts, expressionThreshold)
cts <- cts[highExpr,]

# contrasts to test
contrasts <- data.frame(c('1', '2', '4', '8'), c('0', '0', '0', '0'))
colnames(contrasts) <- c('cond', 'ref')

# construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=~time+data)

# fit all at once
dds$time <- relevel(dds$time, ref='0')
dds <- DESeq(dds)
                                          
apply(contrasts, 1, get_results, dds=dds)
