#! /usr/bin/env Rscript

# Usage: ./plot.R

# figure 3 and associated supplementary figures

library(dplyr)
library(tibble)
library(purrr)

library(ggplot2)
library(RColorBrewer)

library(cowplot)
library(grid)
library(gridExtra)


## input/params

loc <- here::here()
# utils
source(file.path(loc, "paper", "figures", "utils.R", fsep=.Platform$file.sep)) 

# figure 3
figDir <- file.path(loc, "paper", "figures", "3", fsep=.Platform$file.sep) 
suppDir <- file.path(loc, "paper", "figures", "supp", fsep=.Platform$file.sep)

# matched results
rdsDir <- file.path(loc, "paper", "tables", fsep=.Platform$file.sep) 
# main for 0-1-2-4-8
tSet <- '0-1-2-4-8'

# pulseR, GS c('#3182BD', '#31A354')
pal <- c(rev(scales::brewer_pal(pal="Blues")(6))[2],  
         rev(scales::brewer_pal(pal="Greens")(6))[2])


## local functions

corrLabels <- function(n) {
    labels <- c()
    for (c in gsub("\n", " ", n)) {
        components <- unlist(strsplit(c, " "))
        labels <- c(labels, bquote(.(components[1])[.(components[2])]))
    }
    names(labels) <- n
    labels
}


corrLabelsFull <- function(n) {
    labels <- c()
    for (c in gsub("\n", " ", n)) {
        components <- unlist(strsplit(c, " "))
        labels <- c(labels, bquote(.(components[1])[.(components[2])]^.(paste(components[-1:-2], collapse=","))))
    }
    names(labels) <- n
    labels
}


plotCorrd <- function(r, l, n, m=0.6) {
    title <- substitute(paste("Decay rate ", delta, " ", h^-1, " (n=",n, ")"), list(n=n))
    myHeatmap(round(cor(as.matrix(r), use='p'), 2), l, llim=m) + ggtitle(title) +
    theme(plot.title = element_text(hjust = 1.5, vjust=1.5))
}


# functions adapted from https://github.com/dieterich-lab/DesignMetabolicRNAlabeling
plotCIOverlap <- function(cis, title=NULL, leg="none") {
    q <- ggplot(data = cis) +
        geom_linerange(aes(
            ymin = min,
            ymax = max,
            x = i,
            colour = Method
        ), size = 0.2) +
        scale_colour_manual(values = c(pal, "#DCE0E0")) +
        geom_point(aes(
            x = i, y = d),
            colour = "gray65",
            size = 0.05) +
        scale_y_log10(breaks=c(0.1, 1), labels=c(0.10, 1.0)) +
        #xlab("Ordered by rate (from slow to fast)") +
        #ylab("Overlap 95% CI") +
        theme_classic() + 
        theme(legend.position = leg,
              legend.title = element_blank(),
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank(), 
              axis.title.x=element_blank(),
              axis.text.y = element_text(size = 10),
              axis.title.y=element_blank()) +
        guides(color = guide_legend(override.aes = list(size=4, alpha = .9))) +
        ggtitle(title)
    q
}


plotCIOSet <- function(df, title) {
    # use the ERCC decay rate, and determine overlap with its CIs
    df <- df[order(df[,1]),]
    decay <- df[,1]
    l <- list()
    for (data in c("SLAM", "TLS", "TUC")) {
        m <- paste(data, "ERCC", sep = "|")
        d <- df[,grep(m, colnames(df))]
        d <- d[,grep("d$", colnames(d), invert=T)]
        # find overlap
        d$o1 <- pmax(d[[1]], d[[3]])
        d$o2 <- pmin(d[[2]], d[[4]])
        d$o3 <- pmax(d[[1]], d[[5]])
        d$o4 <- pmin(d[[2]], d[[6]])
        # pulseR
        p <- data.frame(
            d = decay,
            min = d$o1,
            max = d$o2,
            i = seq_along(d[,1]),
            Method = "pulseR"
        )
        g <- data.frame(
            d = decay,
            min = d$o3,
            max = d$o4,
            i = seq_along(d[,1]),
            Method = "GS"
        )
        d <- rbind(p, g)
        # pulseR has undefined intervals that are not plotted anyway
        d <- na.omit(d)
        # add extra factor for non-overlapping intervals 
        d$Method <- as.character(d$Method)
        d[(d$min>d$max),]$Method <- "Nonoverlapping"
        d$Method <- factor(d$Method, levels=c("pulseR", "GS", "Nonoverlapping"))
        q <- plotCIOverlap(d, title=data, leg="top")
        l[[data]] <- q
    }
    # combine
    fig <- plot_grid(
            l[["SLAM"]] + theme(legend.position="none"),
            l[["TLS"]] + theme(legend.position="none"),
            l[["TUC"]] + theme(legend.position="none"),
            align = "v",
            ncol = 1)
    legend <- get_legend(l[["SLAM"]] + theme(legend.position = "top",
                                    legend.box.margin = margin(0, 0, 0, 0)) +
                        guides(color = guide_legend(nrow = 1, override.aes = list(size=4, alpha = .9))))
    # create common x and y labels
    y.grob <- textGrob("Overlap 95% CI", 
                       gp=gpar(fontsize=12), rot=90)
    x.grob <- textGrob("Ordered by rate (from slow to fast)", 
                       gp=gpar(fontsize=12))
    t.grob <- NULL
    if (!is.null(title)) {
        title <- paste(paste(unlist(strsplit(title, "-")), collapse=","), "h", sep="")
        t.grob <- textGrob(title, gp=gpar(fontsize=12))
    }
    fig <- plot_grid(legend, fig, ncol = 1, rel_heights = c(.1, 1))
    grid.arrange(arrangeGrob(fig, left = y.grob, bottom = x.grob, top=t.grob))
}


## Call

res <- map(dir(rdsDir, "tbl-0"), ~readRDS(file.path(rdsDir, .x)))
timeSets <- extractTimesFromNames(dir(rdsDir, "tbl-0"))
names(res) <- timeSets

# first plot CI overlap - before removing CIs
supp <- map2(res, names(res), plotCIOSet)

# g <- plot_grid(supp[[1]], supp[[2]], supp[[3]], labels = "auto", ncol = 1)
# ggsave(filename=file.path(suppDir, "figure3_cio.pdf", fsep=.Platform$file.sep),
#        plot = g,
#        width = 180,
#        #height = 210,
#        units = "mm",
#        dpi = 1200)
       
# adjust all names
res <- map(res,
           function(df) {
           d <- df[, grep("d$", colnames(df))]
           n <- gsub(" d$", "", colnames(d))
           n <- gsub("STD", "BSA", n)
           n <- gsub("pulseR ", "pulseR\n", n)
           n <- gsub("GS ", "GS\n", n)
           #n <- unlist(strsplit(n, " "))[seq(1, 2*dim(d)[2], 2)]
           colnames(d) <- n
           d
})

# all together
# names(res) <- NULL
# all_res <- do.call(cbind, res)
# labels <- corrLabelsFull(colnames(all_res))
# heatmap <- plotCorrd(all_res, labels, dim(all_res)[1], m=0.5)
# ggsave(file.path(suppDir, paste("figure3_dall.pdf", sep=""), fsep=.Platform$file.sep), plot=heatmap)

# corr heatmap
supp_heatmaps <- map(res, function(r) {
        labels <- corrLabels(colnames(r))
        plotCorrd(r, labels, dim(r)[1])
    })
names(supp_heatmaps) <- names(res)

# pairs
# we need to write them to disk first, does not work as function...
map2(res, names(res), function(x, y) {
    filen <- file.path(suppDir, paste("pair-", y, ".pdf", sep=""), fsep=.Platform$file.sep)
    pdf(filen, width=16, height=16)
    myPairs(x, 
        lower.panel = lower.panel, 
        upper.panel = myLine,
        cex.labels = 4, cex.axis = 2.5,
        main = "")
    dev.off()    
})

map2(supp_heatmaps, names(supp_heatmaps), function(x, y) {
    filen <- file.path(suppDir, paste("pair-", y, ".pdf", sep=""), fsep=.Platform$file.sep)
    p2 <- ggdraw() + draw_image(filen)
    g <- plot_grid(x, p2, labels = "auto", ncol = 2, nrow = 1, align = 'h')
    # main
    if (y == tSet) {
        save_plot(file.path(figDir, "figure3.pdf", fsep=.Platform$file.sep), g)
    } else {
        save_plot(file.path(suppDir, paste("figure3_d-", y, ".pdf", sep=""), fsep=.Platform$file.sep), g)
    }
})
