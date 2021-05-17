#! /usr/bin/env Rscript


library(topGO)
library(CellPlot)

library(dplyr)
library(purrr)
library(tibble)

library(openxlsx)


# ---------------------------------------------------------

mapping <- "org.Hs.eg.db"
ID <- "Ensembl"
# "classicCount" or "weight01Count" extension class algorithms dealing with the GO graph structure
go.class <- "weight01Count"
# use by default Fisher Test
nodesize <- 10
topNodes <- 250

pvalCutOff <- 0.05 # everywhere 

# ---------------------------------------------------------

## Call

loc <- here::here("dge")
dirloc.out <- file.path(loc, 'results')
contrast <- 'time_8_vs_0.xlsx'
main <- "GO enrichment 8h vs. 0h labeling"

call_cellplot <- function(data) {

    dge <- read.xlsx(file.path(dirloc.out, data, contrast),
                     sheet=2, # use shrunken logFC
                     rowNames=TRUE)
    dge$padj[is.na(dge$padj)] <- 1
    # background - final DGE set
    background <- rownames(dge)
    universe <- background %in% background
    selection <- background %in% rownames(dge[dge$padj<pvalCutOff,])
    relevant.genes <- factor(as.integer(selection[universe]))
    names(relevant.genes) <- background

    topGO.data <- new("topGOdata", 
                      ontology='BP', 
                      allGenes=relevant.genes, 
                      mapping=mapping, 
                      annotationFun=annFUN.org, 
                      nodeSize=nodesize,
                      ID=ID)
    # test statistic
    test.stat <- new(go.class, 
                     testStatistic=GOFisherTest, 
                     name="Fisher")
    # run Fisher test
    sig.groups <- getSigGroups(topGO.data, test.stat)
    # output results
    fisher.results <- GenTable(topGO.data,
                               pvalCutOff=sig.groups,
                               topNodes=topNodes) #length(topGO.data@graph@nodes))
    fisher.results$pvalCutOff <- as.numeric(stringr::str_replace_all(fisher.results$pvalCutOff, "[^0-9e\\-\\.]*", ""))
    fisher.results$LogEnriched <- log2(fisher.results$Significant / fisher.results$Expected)

    ga <- genesInTerm(topGO.data) # GenesAnnotated | list of genes per go-terms
    ga <- ga[fisher.results$GO.ID] # eliminate missing terms
    names(ga) <- NULL
    fisher.results$GenesAnnotated <- ga
    xs <- dge[,c("padj", "log2FoldChange")] # significant stats subset
    xs <- subset(xs, padj < pvalCutOff)
    fisher.results$GenesSignificant <- lapply(fisher.results$GenesAnnotated, intersect, rownames(xs)) # extract genes
    ei.rows <- mclapply(fisher.results$GenesSignificant, function (y) {
    if (length(y)) as.list(xs[y,,drop=FALSE])
    else as.list(rep(NA_real_, length(xs)))
    }, mc.cores = 10)
    ei <- mclapply(names(xs), function(z) {
    lapply(ei.rows, "[[", z)
    }, mc.cores = 10)
    ei <- structure(ei, names = names(xs), row.names = seq(nrow(fisher.results)), class = "data.frame")
    row.names(ei) <- NULL
    fisher.results <- data.frame(fisher.results, ei, stringsAsFactors = FALSE, check.names = FALSE)

    x <- head(fisher.results, 10)
    
    filen <- paste(gsub('.xlsx', '', contrast, fixed = T), "_cellplot.pdf", sep = "")
    pdf(file.path(dirloc.out, data, filen), width=12, height=8)
    cell.plot(x=setNames(x$LogEnrich, x$Term), 
              cells=x$log2FoldChange, 
              main=main, 
              x.mar=c(.5, 0), 
              key.n=7, 
              y.mar=c(.1, 0), 
              cex=1.6, 
              cell.outer=2, 
              bar.scale=.7, 
              space=.2)
    sym.plot(x=setNames(x$LogEnrich, x$Term), 
             cells=x$log2FoldChange, 
             x.annotated = x$Annotated, 
             main=main,
             x.mar=c(.45, 0), 
             key.n=7, 
             cex=1.6, 
             axis.cex=.8, 
             group.cex=.7) 
    dev.off()
}

# lapply(c("slam", "tuc", "tls"), function(data) call_cellplot(data))

call_cellplot('conversion-protocols')

