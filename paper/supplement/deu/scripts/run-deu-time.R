#!/biosw/R/3.6.2/bin/Rscript

## Usage: ./run-deu.R [time]
## time: 1, 2, 4, or 8

# DEXSeq-noaggreg, i.e. without default aggregation of overlapping genes with featureCounts-flat   
#                            
# Isoform prefiltering before construction of counting bins may improve performance...
# but here we just rely on independent filtering implemented in DEXSeq (adapted from DESeq2) - 
# excludes counting bins with expression values that are too low from the testing

library(DEXSeq)

library(dplyr)
library(tibble)
library(purrr)

library(openxlsx)

library("org.Hs.eg.db")

# ---------------------------------------------------------

BPPARAM <- MulticoreParam(38)

# ---------------------------------------------------------

## Functions

getData <- function(data, dataDir, dataFile) {
    dcounts <- read.table(file.path(dataDir, dataFile, fsep=.Platform$file.sep), sep = "\t", header=TRUE, check.names=FALSE)
    dcounts <- dcounts %>% dplyr::select(starts_with("/prj"))
    depth <- length(strsplit(colnames(dcounts)[1], "/", fixed=T)[[1]])
    sampleData <- data.frame(vapply(strsplit(colnames(dcounts), "/", fixed=T), "[", "", depth), stringsAsFactors=FALSE)
    colnames(sampleData) <- "sample"
    sampleData$condition <- vapply(strsplit(vapply(strsplit(sampleData$sample, "_", fixed=T), "[", "", 3), ".", fixed=T), "[", "", 2)
    sampleData$condition <- as.factor(sampleData$condition)
    sampleData$condition <- relevel(sampleData$condition, ref='unlabelled')
    sampleData$time <- as.numeric(gsub("h", "", vapply(strsplit(vapply(strsplit(sampleData$sample, "_", fixed=T), "[", "", 3), ".", fixed=T), "[", "", 1)))
    sampleData$rep <- gsub("^[1-9][0-9][0-9][0-9][0-9][0-9]", "", vapply(strsplit(sampleData$sample, "_", fixed=T), "[", "", 1))
    sampleData$rep <- as.factor(sampleData$rep)
    sampleData$sample <- gsub("[A-Z]$", "", vapply(strsplit(sampleData$sample, "_", fixed=T), "[", "", 1))
    sampleData$sample <- paste(sampleData$sample, "L", sep="")
    sampleData$sample[grepl("unlabelled", sampleData$condition)] <- gsub("L", "U", sampleData$sample[grepl("unlabelled", sampleData$condition)])
    rownames(sampleData) <- sampleData$sample
    sampleData$protocol <- data
    sampleData$protocol <- as.factor(sampleData$protocol)
    list("sampleData"=sampleData, "dcounts"=dcounts)
}

# adapted from https://github.com/vivekbhr/Subread_to_DEXSeq/blob/master/load_SubreadOutput.R
DEXSeqDataSetFromFeatureCounts <- function (countfile, sampleData, 
                                            design = ~sample + exon + condition:exon,
                                            flattenedfile = NULL) 

  {
        # Take a fcount file and convert it to dcounts for dexseq
        message("Adding Exon IDs for DEXSeq")
        countfile %>% dplyr::arrange(Geneid, Start, End) %>% dplyr::select(-c(2:6)) -> dcounts

        colnames(dcounts) <- c("GeneID", rownames(sampleData) )
        id <- as.character(dcounts[,1])
        n <- id
        split(n,id) <- lapply(split(n ,id), seq_along )
        rownames(dcounts) <- sprintf("%s%s%03.f",id,":E",as.numeric(n))
        dcounts <- dcounts[,2:ncol(dcounts)]
        dcounts <- dcounts[substr(rownames(dcounts), 1, 1) != "_", ] #remove _ from beginnning of gene name 
        
        ## get genes and exon names out
        splitted <- strsplit(rownames(dcounts), ":")
        exons <- sapply(splitted, "[[", 2)
        genesrle <- sapply(splitted, "[[", 1)
        
        ## parse the flattened file
        if (!is.null(flattenedfile)) {
                aggregates <- read.delim(flattenedfile, stringsAsFactors = FALSE, 
                                         header = FALSE)
                colnames(aggregates) <- c("chr", "source", "class", "start", 
                                          "end", "ex", "strand", "ex2", "attr")
                aggregates$strand <- gsub("\\.", "*", aggregates$strand)
                aggregates <- aggregates[which(aggregates$class == "exon"), # exonic_part
                                         ]
                aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
                aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", 
                                          aggregates$attr)
                # trim the gene_ids to 255 chars in order to match with featurecounts
                longIDs <- sum(nchar(unique(aggregates$gene_id)) > 255)
                warning(paste0(longIDs, 
                               " aggregate geneIDs were found truncated in featureCounts output"), 
                        call. = FALSE)
                aggregates$gene_id <- substr(aggregates$gene_id,1,255)
                
                transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", 
                                    aggregates$attr)
                transcripts <- strsplit(transcripts, "\\+")
                exonids <- gsub(".*exon_number\\s(\\S+).*", "\\1", # exonic_part_number
                                aggregates$attr)
                exoninfo <- GRanges(as.character(aggregates$chr), IRanges(start = aggregates$start, 
                                                                          end = aggregates$end), strand = aggregates$strand)
                names(exoninfo) <- paste(aggregates$gene_id, exonids, 
                                         sep = ":E")
                
                names(transcripts) <- names(exoninfo) 
                if (!all(rownames(dcounts) %in% names(exoninfo))) {
                        stop("Count files do not correspond to the flattened annotation file")
                }
                matching <- match(rownames(dcounts), names(exoninfo))
                stopifnot(all(names(exoninfo[matching]) == rownames(dcounts)))
                stopifnot(all(names(transcripts[matching]) == rownames(dcounts)))
                dxd <- DEXSeqDataSet(dcounts, sampleData, design, exons, 
                                     genesrle, exoninfo[matching], transcripts[matching])
                return(dxd)
        }
        else {
                dxd <- DEXSeqDataSet(dcounts, sampleData, design, exons, 
                                     genesrle)
                return(dxd)
        }
  
}

# ---------------------------------------------------------

## Call

args <- commandArgs(trailingOnly=TRUE)
data <- as.integer(args[1])

loc <- here::here("deu")
dirloc.out <- file.path(loc, 'results', data)
dir.create(dirloc.out, showWarnings=FALSE)
print(paste("Writing to ", dirloc.out, " ...", sep=""))

dataDir <- file.path(loc, 'tables')
tsv <- "-0h-1h-2h-4h-8h.counts.tsv"

flattenedFile <- file.path(dataDir, 'GRCh38_85_plus_Spikes_Isabel_flat.gtf')

l <- lapply(c("slam", "tuc", "tls"), 
        function(data) getData(data, dataDir, paste(data, tsv, sep="")))
sampleData <- do.call(rbind, lapply(l, `[[`, 1))
dcounts <- do.call(cbind, lapply(l, `[[`, 2))
dcounts <- round(dcounts)

# read exon info 
read.table(file.path(dataDir, paste('slam', tsv, sep=""), fsep=.Platform$file.sep), skip=1, header=TRUE) %>% 
    dplyr::select(c(1:6)) -> exoninfo
countfile <- cbind(exoninfo, dcounts)
    
# use all samples from split files, i.e. labelled and unlabelled for the conversion protocols
designFullModel    =  ~ sample + exon + protocol:exon + condition:exon
designReducedModel =  ~ sample + exon + protocol:exon 


dxd <- DEXSeqDataSetFromFeatureCounts(countfile, 
                                      sampleData, 
                                      design=designFullModel, 
                                      flattenedfile=flattenedFile) 
# select time point
dxd <- dxd[,dxd$time==data]

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, 
                           formula=designFullModel,
                           BPPARAM=BPPARAM)
dxd <- testForDEU(dxd, 
                  reducedModel=designReducedModel, 
                  fullModel=designFullModel,
                  BPPARAM=BPPARAM)
dxd <- estimateExonFoldChanges(dxd, 
                               fitExpToVar="condition", 
                               denominator='unlabelled',
                               BPPARAM=BPPARAM)

filen <- paste0(data, "_dxd.rds", sep="")
filen <- file.path(dirloc.out, filen, fsep=.Platform$file.sep)
saveRDS(dxd, file=filen)
                               
dxr <- DEXSeqResults(dxd)
dxr <- dxr %>% as.data.frame() %>%
            mutate(transcripts = sapply(transcripts, toString))
            
dxr$symbol <- mapIds(org.Hs.eg.db,
                     keys=dxr$groupID,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
dxr <- dxr[,c(ncol(dxr),1:(ncol(dxr)-1))]                        
                         
# write to disk 
wb <- createWorkbook()
addWorksheet(wb, sheetName='LvsU')
writeDataTable(wb, sheet=1, x=dxr)
filen <- paste0(data, "_DEU_LRT.xlsx", sep="")
filen <- file.path(dirloc.out, filen, fsep=.Platform$file.sep)
saveWorkbook(wb, filen, overwrite = TRUE)
