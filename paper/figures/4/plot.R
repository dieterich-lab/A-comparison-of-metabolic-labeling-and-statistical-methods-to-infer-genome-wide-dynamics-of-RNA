#!/biosw/R/3.5.1/bin/Rscript

# Usage: ./plot.R

# figure 4 and associated supplementary figures

library(dplyr)
library(tibble)
library(purrr)

library(UpSetR)

library(ggplot2)
library(RColorBrewer)

library(cowplot)
library(grid)
library(gridExtra)

## input/params

loc <- here::here()
# utils
source(file.path(loc, "paper", "figures", "utils.R", fsep=.Platform$file.sep)) 

# figure 4
figDir <- file.path(loc, "paper", "figures", "4", fsep=.Platform$file.sep) 
suppDir <- file.path(loc, "paper", "figures", "supp", fsep=.Platform$file.sep)

# matched results
rdsDir <- file.path(loc, "paper", "tables", fsep=.Platform$file.sep) 
# main for 0-1-2-4-8
tSet <- '0-1-2-4-8'

# we need all the results and the counts for each method/protocol
# pulseR
pulseDir <- file.path(loc, "pulseRTc", "results", fsep=.Platform$file.sep) 
dataSets <- c("ercc", "std", "slam", "tls", "tuc")
pulseFit <- "pulsefit"
pulseCis <- "pulsecis"
# GRAND-SLAM
gsDir <- file.path(loc, "grand-slam", fsep=.Platform$file.sep)
gsFit <- "all"

# pulseR, GS, pulseR, GS - box plot order
pal <- c('#3182BD', '#31A354', '#6BAED6', '#74C476')
col_order <- c('pulseR shared', 'GS shared', 'pulseR only', 'GS only')

# Methods "#000000", "#999999", # ERCC, STD
# slam, tls, tuc, ercc, std
mpal <- c('#abc3cc', '#80a2ab', '#658087', "#000000", "#999999")

# bottom to top - upset
setsm <- c("ERCC pulseR", "BSA pulseR",
          "TUC GS", "TUC pulseR",
          "TLS GS", "TLS pulseR",
          "SLAM GS", "SLAM pulseR")
          
metadata <- data.frame(sets=setsm, 
                       hue=c("ercc", "bsa",
                             "gs", "pr", 
                             "gs", "pr",
                             "gs", "pr"))
                             
                             
## local functions
                
myLine <- function(x, y, ...){
    par(mar=c(5, 5, 4, 2))
    smoothScatter(x, y, 
                  colramp = shadesOfGrey,
                  cex = 3, cex.axis = 2,
                  xaxt='n', yaxt = 'n',
                  ...)
    axis(1, at = c(0, 2, 4), labels = 10^c(0, 2, 4), cex.axis = 2)
    axis(2, at = c(-1, 0, 1), labels = 10^c(-1, 0, 1), cex.axis = 2)
    log10.axis(1, at=c(0, 2, 4))
    log10.axis(2, at=c(-1, 0, 1))
    abline(a = 0,b = 1, ...)
}


plotMyLine <- function(x, y, m) {
    myLine(x, y, xlab = expression("GS (log"[10]*" read count)"), 
            ylab = expression("pulseR (log"[10]*" read count)"), 
            cex.lab=2)
    title(m, cex.main=2.5, adj = 0)
}


# redefined from package UpSetR
fromList_ <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
      x <- as.vector(match(elements, x))
      }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names
  row.names(data) <- elements
  return(data)
  }


# supp. figure
getPlot <- function(df, main, mycolor1, mycolor2){
    yterms <- unique(unlist(lapply(df, function(x){
      x <- subset(x, pvalCutOff <= 0.05) # already filtered anyway...!
      x <- x[order(x$pvalCutOff),] # order by pvalCutOff, but plot LogEnriched
      head(x, 5)$GO.ID
    })))
    par(mar = c(1,1,1,10)) # all hard coded
    myHistogram(df, min.genes = 1, max.genes = 1e10, go.selection = yterms, 
             main = main, 
             axis.cex = 1.5, lab.cex = 2, main.cex = 2.5, 
             mycolor1 = mycolor1, mycolor2 = mycolor2, fixed.lim=8) # hard coded fixed limit
}


## call

# load pulseR estimates
pulseFiles <- lapply(tSet, function(x) paste(paste(pulseFit, paste(x, collapse = "-"), sep = "-"), "rds", sep = "."))
pulseFiles <- file.path(pulseDir , rep(dataSets, each = length(pulseFiles)), pulseFiles, fsep=.Platform$file.sep)

pulseList <- map(pulseFiles, function(.id) {
  data <- toupper(rev(unlist(strsplit(.id, "/")))[2])
  pulse <- readRDS(.id)
  pulsedf <- as.data.frame(pulse$fit$d)
  cis <- readRDS(gsub(pulseFit, pulseCis, .id))
  cis <- as.data.frame(cis)
  df <- cbind(pulsedf, cis)
  colnames(df) <- c("d", "d.min", "d.max")
  df$gene <- rownames(pulse$pd$counts)
  # adjust name
  if (data == 'STD') { data <- 'BSA'}
  df$Method <- data 
  df$Workflow <- 'pulseR'
  df
})

pulseNames <- map(pulseFiles, function(.id) {
  toupper(rev(unlist(strsplit(.id, "/")))[2])
})
names(pulseList) <- pulseNames

# load GRAND-SLAM estimates (from get_estimates.R) - includes all for selected time set 
gsTables <- file.path(gsDir, "tables", fsep=.Platform$file.sep)
gsNames <- map(dir(gsTables, paste("[slam|tls|tuc]", gsub('0', '', tSet), ".rds", sep="")), ~toupper(unlist(strsplit(.x, "-"))[1]))
gsList <- map(dir(gsTables, paste("[slam|tls|tuc]", gsub('0', '', tSet), ".rds", sep="")), ~readRDS(file.path(gsTables, .x)))
names(gsList) <- gsNames

gsList <- map2(gsList, names(gsList), function(df, n) {
  colnames(df) <- c("d", "d.min", "d.max")
  # remove where we have no predictions - GS include NAs! - only for d
  df <- df[!is.na(df$d),]
  df$gene <- rownames(df)
  rownames(df) <- NULL
  df$Method <- n
  df$Workflow <- 'GS'
  df
})

# combine into tidy format, add summary statistics per gene
tidy <- do.call(rbind, c(pulseList, gsList))
s <- tidy %>%
       group_by(gene) %>%
       dplyr::summarize(d.median = median(d), d.mean = mean(d), n = n())
tidy <- merge(tidy, s, by='gene', all=TRUE)
tidy$w <- (tidy$d.max - tidy$d.min) # use absolute difference, i.e. actual width
s <- tidy %>%
       group_by(gene) %>%
       dplyr::summarize(w.median = median(w), w.mean = mean(w))
tidy <- merge(tidy, s, by='gene', all=TRUE)

# first create the upset matrix/table
select <- tidyr::expand_grid(unique(tidy$Method),unique(tidy$Workflow))
select$set <- as.character(interaction(select,sep=" "))
select <- select[!select$set %in% c("BSA GS", "ERCC GS"),]

sets <- list()
for (s in unique(select$set)) {
    sets[[s]] <- tidy[tidy$Method==unlist(strsplit(s, " "))[1] & tidy$Workflow==unlist(strsplit(s, " "))[2],]$gene
}
upsetdf <-  fromList_(sets)
upsetdf$gene <- rownames(upsetdf)
upsetdf <- merge(upsetdf, unique(tidy[,c("gene", "d.median", "w.median")]), by='gene', all.x=T, all.y=F)

# upset - save figure and reload
assignInNamespace("BoxPlotsPlot", myBoxPlotsPlot, ns="UpSetR")
assignInNamespace("BaseBoxPlot", myBaseBoxPlot, ns="UpSetR")

# find cardinality, name sets, and remove intersections with less than 10 for visualization
setNames <- names(upsetdf)[2:(dim(upsetdf)[2]-2)]
n <- length(setNames)
z <- rep(0,n)
# most are empty...
combinations <- t(apply(combn(0:n,1), 2, function(.id) {z[.id]=1; z}))
for (i in 2:n) {
    combinations <- rbind(combinations, t(apply(combn(0:n,i), 2, function(.id) {z[.id]=1; z})))
}
combinations <- combinations[!duplicated(combinations),]

selectedCols <- upsetdf[,setNames]
sets <- apply(combinations, 1, function(.id) {
    selectedRows <- selectedCols[apply(selectedCols, 1, 
                                   function(r) {all(r==as.integer(.id))}), ]
    upsetdf[rownames(upsetdf) %in% rownames(selectedRows),]$gene
})
names(sets) <- apply(combinations, 1, function(.id) {
    paste(setNames[as.logical(.id)], collapse="+")
})
sets <- sets[lapply(sets,length)>0]
small <- sets[lapply(sets,length)<=10]
upsetdfLarge <- upsetdf
for (s in small) {
    upsetdfLarge <- upsetdfLarge[!upsetdfLarge$gene %in% s,] 
}

filen <- file.path(suppDir, paste("upset", tSet, "pdf", sep="."), fsep=.Platform$file.sep)
pdf(filen, width=20, height=16, onefile=FALSE) # https://github.com/hms-dbmi/UpSetR/issues/115
upset(upsetdfLarge,      
      sets = setsm, 
      order.by = "degree", keep.order = TRUE, 
      mainbar.y.label = "Genes in common",
      sets.x.label = "Total genes",
      sets.bar.color=c(mpal[4], mpal[5],
                       mpal[3], mpal[3],
                       mpal[2], mpal[2],
                       mpal[1], mpal[1]),
      set.metadata = list(data = metadata, plots = list(list(type = "matrix_rows", 
      column = "hue", 
      colors = c(gs = pal[2], pr = pal[1], 
                 ercc = mpal[4], bsa = mpal[5]),
                 alpha = 0.6))),
      text.scale = 3, point.size = 5, line.size = 1.25, mb.ratio = c(0.55,0.45),
      boxplot.summary = c("d.median", "w.median")
)
dev.off()

# box plot

# find "sets" containing genes that are shared between (resp. unique to) pulseR and GS, but wrt each protocol only
# keep these for the supp figures
goSets <- list()
tidy$Set <- NA
# pulseR or GS shared - SLAM
g <- upsetdf[upsetdf[['SLAM pulseR']]==1 & upsetdf[['SLAM GS']]==1,]$gene
tidy$Set[tidy$gene %in% g & tidy$Method == "SLAM" & tidy$Workflow == 'pulseR'] <- 'pulseR shared' 
tidy$Set[tidy$gene %in% g & tidy$Method == "SLAM" & tidy$Workflow == 'GS'] <- 'GS shared' 
# pulseR only
g <- upsetdf[upsetdf[['SLAM pulseR']]==1 & upsetdf[['SLAM GS']]==0,]$gene
tidy$Set[tidy$gene %in% g & tidy$Method == "SLAM"] <- 'pulseR only' 
goSets[['SLAM pulseR']] <- tidy[tidy$Set == 'pulseR only' & tidy$Method == "SLAM",]$gene
# GS only
g <- upsetdf[upsetdf[['SLAM pulseR']]==0 & upsetdf[['SLAM GS']]==1,]$gene
tidy$Set[tidy$gene %in% g & tidy$Method == "SLAM"] <- 'GS only'
goSets[['SLAM GS']] <- tidy[tidy$Set == 'GS only' & tidy$Method == "SLAM",]$gene

# pulseR or GS shared - TLS
g <- upsetdf[upsetdf[['TLS pulseR']]==1 & upsetdf[['TLS GS']]==1,]$gene
tidy$Set[tidy$gene %in% g & tidy$Method == "TLS" & tidy$Workflow == 'pulseR'] <- 'pulseR shared' 
tidy$Set[tidy$gene %in% g & tidy$Method == "TLS" & tidy$Workflow == 'GS'] <- 'GS shared' 
# pulseR only
g <- upsetdf[upsetdf[['TLS pulseR']]==1 & upsetdf[['TLS GS']]==0,]$gene
tidy$Set[tidy$gene %in% g & tidy$Method == "TLS"] <- 'pulseR only' 
goSets[['TLS pulseR']] <- tidy[tidy$Set == 'pulseR only' & tidy$Method == "TLS",]$gene
# GS only
g <- upsetdf[upsetdf[['TLS pulseR']]==0 & upsetdf[['TLS GS']]==1,]$gene
tidy$Set[tidy$gene %in% g & tidy$Method == "TLS"] <- 'GS only' 
goSets[['TLS GS']] <- tidy[tidy$Set == 'GS only' & tidy$Method == "TLS",]$gene

# pulseR or GS shared - TUC
g <- upsetdf[upsetdf[['TUC pulseR']]==1 & upsetdf[['TUC GS']]==1,]$gene
tidy$Set[tidy$gene %in% g & tidy$Method == "TUC" & tidy$Workflow == 'pulseR'] <- 'pulseR shared' 
tidy$Set[tidy$gene %in% g & tidy$Method == "TUC" & tidy$Workflow == 'GS'] <- 'GS shared' 
# pulseR only
g <- upsetdf[upsetdf[['TUC pulseR']]==1 & upsetdf[['TUC GS']]==0,]$gene
tidy$Set[tidy$gene %in% g & tidy$Method == "TUC"] <- 'pulseR only' 
goSets[['TUC pulseR']] <- tidy[tidy$Set == 'pulseR only' & tidy$Method == "TUC",]$gene
# GS only
g <- upsetdf[upsetdf[['TUC pulseR']]==0 & upsetdf[['TUC GS']]==1,]$gene
tidy$Set[tidy$gene %in% g & tidy$Method == "TUC"] <- 'GS only' 
goSets[['TUC GS']] <- tidy[tidy$Set == 'GS only' & tidy$Method == "TUC",]$gene

# remove BSA/ERCC
tidy <- tidy[!tidy$Method %in% c("ERCC", "BSA"),]

# redorder as factor 
col_order <- c('pulseR shared', 'GS shared', 'pulseR only', 'GS only')
tidy$Set <- factor(tidy$Set, levels = col_order)

q <- ggplot(tidy, aes(x=Set, y=d, fill=Set)) +
    #geom_jitter(aes(color=Set), alpha = 0.3, width = .25) +
    geom_boxplot(alpha = 0, colour = "black", show.legend=FALSE, fatten = .5, lwd = .5) +
    facet_wrap( ~ Method) +
    scale_colour_manual(values = pal) +
    ylab(expression(paste(delta, ", ", h^-1))) +
    stat_summary(
        fun.data = statBoxData, 
        geom = "text", 
        size = 2,
        position = position_dodge(width = 0.75)
    ) +
    ggpubr::stat_compare_means(method="wilcox.test", 
                            comparisons = list(c("pulseR shared", "GS shared"), c("pulseR only", "GS only"),
                            c("pulseR shared", "pulseR only"), c("GS shared", "GS only")), 
                            label = "p.signif",
                            paired=FALSE,
                            size = 2) +
    theme_classic() + 
            theme(legend.position = "bottom",
                legend.title = element_blank(),
                legend.text=element_text(size=10),
                axis.ticks.x = element_blank(), 
                axis.text.x = element_blank(), 
                axis.title.x=element_blank(),
                axis.text.y = element_text(size = 10), 
                axis.title.y = element_text(size = 10),
                strip.text = element_text(size = 10),
                strip.background = element_blank()) + #panel.spacing = unit(3, "lines")
    guides(color = guide_legend(override.aes = list(size=2, alpha = .9))) +
    scale_x_discrete(expand=c(0.4,0))
# for jitter points only, downsample to reduce figure size
small <- tidy[tidy$Set %in% c("pulseR only", "GS only"),]
large <- sample_frac(tidy[tidy$Set %in% c("pulseR shared", "GS shared"),], 0.25, weight = d)
q <- q + geom_jitter(rbind(small, large), mapping = aes(x=Set, y=d, color=Set), alpha = 0.025, width = .3, size= 0.5)

# scatter plots of raw counts for selected subsets

gsOnly <- map2(sets, names(sets), function(s, n) {
    if (grepl('pulseR', n)) { }
    else { s }
})
gsOnly <- gsOnly[lapply(gsOnly,length)>0]

# pulseR results - median (across protocols) of mean total counts (across samples, per protocol)
fcList <- map(dataSets[grep("ercc", dataSets, invert=T)], ~readRDS(file.path(pulseDir, .x, "data", "counts.rds")))
names(fcList) <- toupper(dataSets[grep("ercc", dataSets, invert=T)])
fcList <- map2(fcList, names(fcList), function(cts, data) {
    if (data == 'STD') {
        # we need to total samples... ands remove ERCC
        cts <- cts[!grepl("ERCC-", rownames(cts)),]
        conditions <- readRDS(file.path(pulseDir, "std", "data", "samples.rds"))
        cts <- cts[,colnames(cts) %in% conditions[conditions$fraction=='total',]$sample]
    } else { # slam, tls, or tuc
        conditions <- data.frame(colnames(cts))
        names(conditions) <- "sample"
        conditions$map <- gsub('L|U', '', conditions$sample)
        totals <- map(unique(conditions$map), function(.id) {
        i <- conditions$sample[conditions$map == .id]
        apply(cts[,i], 1, sum)
        })
        cts <- as.data.frame(do.call(cbind, totals))
        colnames(cts) <- unique(conditions$map)
    }
    # add mean gene count
    cts$mean <- rowMeans(cts)
    cts
})

allMeans <- map(fcList, function(.id) { .id[,c("mean"),drop=FALSE] } )
allMeans <- do.call("cbind", allMeans)
allMeans$median <- robustbase::rowMedians(as.matrix(allMeans)) 
fcList[['ALL']] <- allMeans

# GRAND-SLAM
gsFiles <- file.path(gsDir , dataSets[grep("slam|tuc|tls", dataSets)], gsFit, paste(gsFit, "tsv", sep="."), fsep=.Platform$file.sep)
gsNames <- map(gsFiles, ~toupper(rev(unlist(strsplit(.x, "/")))[3]))
gsList <- map(gsFiles, function(.id) {
    # swe only need to read count...
    # XXX Readcount 	The total number of reads mapped to this gene in condition XXX
    gs <- read.table(.id, sep='\t', header=T, row.names = 1)
    gs <- gs[,grepl('Readcount', colnames(gs))]
    gs$mean <- rowMeans(gs)
    gs
})
names(gsList) <- gsNames

allMeans <- map2(gsList, names(gsList), function(.id, n) {
    .id[,c("mean"),drop=FALSE] %>%
        dplyr::rename(!!n:=mean) %>%
        rownames_to_column(var = "gene")
} )
allMeans <- plyr::join_all(allMeans, by = 'gene', type = 'full') %>%
                 column_to_rownames("gene")
allMeans$median <- robustbase::rowMedians(as.matrix(allMeans), na.rm=TRUE) 
gsList[['ALL']] <- allMeans

# cannot use as function... write each plot to disk first
plot_figs <- function(sets, x, y, basename) {
  qs <- map2(sets, names(sets), function(.x, .y) {
      filen <- file.path(suppDir, paste(basename, "scatter", .y, tSet, "png", sep="."), fsep=.Platform$file.sep)
      png(filen)
      plotMyLine(log10(x[[.y]][rownames(x[[.y]]) %in% .x,]$mean),
                            log10(y[[.y]][rownames(y[[.y]]) %in% .x,]$mean), .y)
      dev.off()
  })
}

# compare the GS-specific gene sets
res <- gsOnly[c('SLAM GS', 'TLS GS', 'TUC GS')]
names(res) <- c('SLAM', 'TLS', 'TUC')
plot_figs(res, gsList, fcList, "gs-excl")
pngfiles <- list.files(path = suppDir,
                       recursive = TRUE,
                       pattern = "^gs-excl.scatter.*png$",
                       full.names = T)
                       
gsAll <- unlist(gsOnly)
filen <- file.path(suppDir, paste("scatter-gs-all", tSet, "png", sep="."), fsep=.Platform$file.sep)
png(filen) #, width = 8, height = 8, units = 'in', res = 1200)
plotMyLine(log10(gsList[['ALL']][rownames(gsList[['ALL']]) %in% gsAll,]$median),
                            log10(fcList[['ALL']][rownames(fcList[['ALL']]) %in% gsAll,]$median), 'All GS intersections')
dev.off()                    

# 4 plot panel
pngfiles <- c(pngfiles, filen)
magick::image_read(pngfiles) %>%
    magick::image_montage(tile = "2x2", geometry = "x500") %>%
    magick::image_convert("pdf") %>%
    magick::image_write(format = "pdf", 
                        path = file.path(suppDir, paste("scatter-panel-final", tSet, "pdf", sep="."), fsep=.Platform$file.sep),
                        quality = 1200)

# final panel
filen <- file.path(suppDir, paste("upset", tSet, "pdf", sep="."), fsep=.Platform$file.sep)
p1 <- ggdraw() + draw_image(filen)

filen <- file.path(suppDir, paste("scatter-panel-final", tSet, "pdf", sep="."), fsep=.Platform$file.sep)
p2 <- ggdraw() + draw_image(filen)

title <- ggdraw() + draw_label("Genes used by\nGRAND-SLAM only")
p22 <- plot_grid(title, p2, ncol=1, rel_heights=c(0.25, 1)) # rel_heights values control title margins

g1 <- plot_grid(p1, p22, labels = c("a", "b"), rel_widths = c(1, .55), ncol = 2)
g <- plot_grid(g1, q, labels = c("", "c"), rel_heights = c(1, .55), nrow = 2)
ggsave(filename=file.path(figDir, "figure4.pdf", fsep=.Platform$file.sep),
       plot = g,
       width = 180,
       #height = 210,
       units = "mm",
       dpi = 1200)

# supp - enrichment on selected gene sets       
       
library('org.Hs.eg.db')

# use topGO for gene set enrichment
library(topGO)

library(wordcloud)
       
# using annFUN.org
ontologies.list <- c('BP', 'CC', 'MF')
mapping <- "org.Hs.eg.db"
ID <- "Ensembl"

# "classicCount" or "weight01Count" extension class algorithms dealing with the GO graph structure
go.class <- "weight01Count"
# use by default Fisher Test
nodesize <- 10
topNodes <- 250

pvalCutOff <- 0.05

# only biological process
ont <- 'BP'
       
# to define the background, we use the featureCounts table
# and filter based on median-mean count
background <- rownames(fcList[['ALL']][fcList[['ALL']]$median > 50,])
universe <- background %in% background

# histogram for SLAM, TLS, and TUC GS vs pulseR 2 colors
res <- map2(goSets, names(goSets), function(.x, .y) {
            selection <- background %in% .x
            relevant.genes <- factor(as.integer(selection[universe]))
            names(relevant.genes) <- background
                
            topGO.data <- new("topGOdata", 
                            ontology=ont, 
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
            fisher.results$Cluster <- .y
            fisher.results
        })
names(res) <- names(goSets) 

# we need to "tweak" the figure...

ids <- c('GO:0002313',
         'GO:0035235',
         'GO:0030512',
         'GO:0007156',
         'GO:0006620',
         'GO:1902902',
         'GO:0010984')
terms <- c('mature B cell diff. involved in IR',
           'iGluR signaling pathway',
           'neg. reg. of TGFB receptor sign. path.',
           'homophilic cell adhesion via CAMs',
           'posttr. prot. targeting to ERM',
           'neg. reg. of autophagosome assembly',
           'reg. of lipoprot. particle clearance')
# adjust manually top names for pretty display
res <- lapply(res, function(.x) {
    for (j in 1:length(ids)) {
        if (ids[j] %in% unique(.x$GO.ID)) {
            .x[.x$GO.ID==ids[j],]$Term <- terms[j]
        }
    }
    .x
})
# arbitrary max.
res[['TUC GS']][res[['TUC GS']]$GO.ID=='GO:0034454',]$LogEnriched <- 8
res[['TUC GS']][res[['TUC GS']]$GO.ID=='GO:0001946',]$LogEnriched <- 8

# define the order as before pulseR then GS
res <- res[c('SLAM pulseR', 'TLS pulseR', 'TUC pulseR', 'SLAM GS', 'TLS GS', 'TUC GS')]

filen <- file.path(suppDir, paste("enrichment", tSet, "pdf", sep="."), fsep=.Platform$file.sep)
pdf(filen, width=14, height=7)
getPlot(res, "GO enrichment\nbiological process", pal[1], pal[2])
dev.off()


# word cloud of biotypes for exclusive pulseR and GS sets across all methods

# redefine for convenience
pal <- list('pulseR'=rev(scales::brewer_pal(pal="Blues")(6))[2],
         'GS'=rev(scales::brewer_pal(pal="Greens")(6))[2])
         
pulseROnly <- map2(sets, names(sets), function(s, n) {
    if (grepl('GS', n)) { }
    else { s }
})
pulseROnly <- pulseROnly[lapply(pulseROnly,length)>0]
pulseRAll <- unlist(pulseROnly)
geneSets <- list(gsAll, pulseRAll)
names(geneSets) <- c('GS exclusive', 'pulseR exclusive')

# annotation file 

# annotations
aloc <- file.path(loc, "pulseRTc", "workflow", "tables", fsep=.Platform$file.sep) 
ann <- read.csv(file.path(aloc, 'featureCount.bed.gz', fsep=.Platform$file.sep), 
                sep="\t",
                check.names=FALSE,
                stringsAsFactors=FALSE)
names(ann) <- stringr::str_replace_all(colnames(ann), "#", "")

map2(geneSets, names(geneSets), function(.x, .y) {
    name <- unlist(strsplit(.y, " "))[1]
    filen <- file.path(suppDir, paste("cloud", name, "pdf", sep="."), fsep=.Platform$file.sep)
    set <- ann[ann$id %in% .x,]
    freq <- as.data.frame(table(set$biotype))
    freq$Var1 <- gsub("_", " ", freq$Var1)
    set.seed(1234)
    pdf(filen, width=5, height=5)
    par(mar=c(0, 0, 0, 0))
    wordcloud(words = freq$Var1, 
              freq = freq$Freq, 
              min.freq = 2,
              random.order = FALSE, 
              rot.per = 0,
              scale=c(2,.5),
              colors=pal[[name]])
    dev.off()
})

# put together 
filen <- file.path(suppDir, paste("enrichment", tSet, "pdf", sep="."), fsep=.Platform$file.sep)
p1 <- ggdraw() + draw_image(filen)

filen <- file.path(suppDir, paste("cloud.pulseR.pdf", sep="."), fsep=.Platform$file.sep)
p2 <- ggdraw() + draw_image(filen)

filen <- file.path(suppDir, paste("cloud.GS.pdf", sep="."), fsep=.Platform$file.sep)
p3 <- ggdraw() + draw_image(filen)

# this doesn't do...
# we still need to rescale the final figure...
g2 <- plot_grid(p2, p3, ncol = 2)
g <- plot_grid(p1, g2, labels = "auto", rel_heights = c(1, .45), nrow = 2)
ggsave(filename=file.path(suppDir, "figureSupp4.pdf", fsep=.Platform$file.sep),
       plot = g,
       width = 180,
       #height = 210,
       units = "mm",
       dpi = 1200)
