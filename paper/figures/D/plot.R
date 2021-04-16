#!/biosw/R/3.5.1/bin/Rscript

# Usage: ./plot.R

# figures associated with supplementary file "Estimating variance and optimal labeling time using the pulseR workflow."

# plot diagonal term of the FIM (standard design and RTC, incl. overdispersion)
# computed at estimated parameters

library(dplyr)
library(tibble)
library(purrr)
library(stringr)

library(ggplot2)
library(RColorBrewer)

library(cowplot)
library(grid)
library(gridExtra)

## Adapted from https://github.com/dieterich-lab/DesignMetabolicRNAlabeling

## input/params

loc <- here::here()
# utils
source(file.path(loc, "paper", "figures", "utils.R", fsep=.Platform$file.sep)) 

# figure 5
figDir <- file.path(loc, "paper", "figures", "D", fsep=.Platform$file.sep) 

# pulseR results only 
pulseDir <- file.path(loc, "pulseRTc", "results", fsep=.Platform$file.sep) 

dataSets <- c("ercc", "std", "slam", "tls", "tuc")
modelStr <- "pulse"
basename <- "figure5"

# main for 0-1-2-4-8
tSet <- '0-1-2-4-8'

# labeled - unlabeled
pal <- c("#999999", "#f0f0f0")
# slam, tls, tuc, ercc, std
mpal <- list('slam'='#abc3cc', 'tls'='#80a2ab', 'tuc'='#658087', 'ercc'="#000000", 'std'="#999999")

## local functions

getIFromDataBSA <- function(fit, normCoef, t, condition) {
    if (condition == "total") return (NULL)
    fun <-c("pull_down" = I_labelledBSA, "flow_through" = I_unlabelledBSA)[condition]
    res <- data.frame(
        I = fun[[1]](list(
        mu = exp(fit$mu)*normCoef,
        d = fit$d,
        k = fit$size,
        t = t
        )),
        t = t,
        condition = c("pull_down" = "Eluate",
        "flow_through" = "Supernatant")[condition],
        stringsAsFactors = FALSE
        )
    res$id <- rownames(res)
    res$d <- fit$d
    res$mu <- fit$mu
    res
}


getIFromDataRTC <- function(fit, normCoef, t, condition) {
    if (t == 0) return (NULL)
    fun <-c("labelled" = I_labelledRTC, "unlabelled" = I_unlabelledRTC)[condition]
    res <- data.frame(
        I = fun[[1]](list(
        mu1 = (exp(fit$mu1) + exp(fit$mu2) + exp(fit$mu3))*normCoef,
        mu2 = exp(fit$mu2)*normCoef,
        mu3 = exp(fit$mu3)*normCoef,
        d = fit$d,
        k = fit$size,
        t = t
        )),
        t = t,
        condition = c("labelled" = "Labeled",
        "unlabelled" = "Unlabeled")[condition],
        stringsAsFactors = FALSE
        )
    res$id <- rownames(res)
    res$d <- fit$d
    res
}


plotFIForGenesBSA <- function(o, r) {
    q <- ggplot(data = sample_frac(o, 0.5, weight = I)) +
        geom_point(aes(x = d, y = d^2*I, colour = condition), stroke = 0, size = .18) +
        facet_wrap(~ t, scales = "free_y", labeller = as_labeller(function(x) paste(x, " h"))) +
        ylab(expression(paste(I[paste(delta, delta)] %.% delta^2))) +
        xlab(expression(paste(delta, ", ", h^-1))) +
        scale_x_log10() +
        scale_y_continuous(trans = "sqrt") +
        #scale_color_brewer(palette = pal, name = "") +
        scale_colour_grey(, name = "") +
        guides(colour = guide_legend(override.aes = aes(size = 1))) +
        annotation_logticks(sides = 'b') +
        geom_line(data = o %>% filter(d < 1),
            aes(x = d, y = t^2*r$fit$size*d^2), linetype = "42") +
        geom_line(data = o %>% filter(d < 1),
            aes(x = d, y = t^2*r$fit$size*d^2 * exp(-2*d*t) / (1 - exp(-d*t))^2), linetype = "42") +
        theme_classic() + 
        theme(strip.background = element_blank(),
              strip.placement = "outside",
              strip.text.x = element_text(size = 12),
              legend.position = c(.73,1.2),
              legend.justification = c(1,1),
              axis.title.y = element_text(family = "serif"),
              text = element_text(size = 12),
              axis.text = element_text(size = 12)) +
        guides(colour = guide_legend(override.aes = list(size=2)))
    q
}


# we plot the asymptotic curve for "eluate" (std)
# for RTC, this does not depend on mu1, in the limit of mu1 >> mu2 + mu3
# and the unlabelled FIM diagonal term goes to zero
plotFIForGenesRTC <- function(o, r) {
    q <- ggplot(data = sample_frac(o, 0.4, weight = I)) +
        geom_point(aes(x = d, y = d^2*I, colour = condition), stroke = 0, size = .18) +
        facet_wrap(~ t, scales = "free_y", labeller = as_labeller(function(x) paste(x, " h"))) +
        ylab(expression(paste(I[paste(delta, delta)] %.% delta^2))) +
        xlab(expression(paste(delta, ", ", h^-1))) +
        scale_x_log10() +
        scale_y_continuous(trans = "sqrt") +
        #scale_color_brewer(palette = pal, name = "") +
        scale_colour_grey(, name = "") +
        guides(colour = guide_legend(override.aes = aes(size = 1))) +
        annotation_logticks(sides = 'b') +
        geom_line(data = o %>% filter(d < 1),
                  aes(x = d, y = t^2*r$fit$size*d^2 * exp(-2*d*t) / (1 - exp(-d*t))^2), linetype = "42") +
        theme_classic() + 
        theme(strip.background = element_blank(),
              strip.placement = "outside",
              strip.text.x = element_text(size = 12),
              legend.position = c(.175,1),
              legend.justification = c(1,1),
              axis.title.y = element_text(family = "serif"),
              text = element_text(size = 12),
              axis.text = element_text(size = 12)) +
        guides(colour = guide_legend(override.aes = list(size=2)))
    q
}


fiForGenesBSA <- function(r) {
    norms <- pulseR:::getNorms(r$pd, r$fit$normFactors)
    norms <- apply(norms, 2, max)
    o <- pmap(list(normCoef = as.list(norms),
                   t = as.list(r$pd$conditions$time),
                   condition = as.character(r$pd$conditions$fraction)),
              getIFromDataBSA,
              fit = r$fit)
    o <- bind_rows(o)
    plotFIForGenesBSA(o, r)
}


fiForGenesRTC<- function(r) {
    norms <- r$pd$depthNormalisation
    o <- pmap(list(normCoef = as.list(norms),
                   t = as.list(r$pd$conditions$time),
                   condition = as.character(r$pd$conditions$fraction)),
              getIFromDataRTC,
              fit = r$fit)
    o <- bind_rows(o)
    plotFIForGenesRTC(o, r)
}


plotCIOrdered <- function(fit, cis, col) {
  o <- order(fit$d)
  cis <- naToBorders(cis, c(1e-3, 2))
  d <- data.frame(
    d = fit$d[o],
    min=cis[o,1],
    max=cis[o,2],
    i = seq_along(cis[,1])
  )
  ggplot(data = d) +
    geom_linerange(aes(
      ymin = min,
      ymax = max,
      x = i
    ),
    colour = col, alpha = .25, size = .25) +
    geom_point(aes(
      x = i, y = d),
      colour = "gray65", size = .25) +
    scale_y_log10() +
    xlab("Ordered by rate (from slow to fast)") +
    ylab(expression(paste(delta, ",  ", h^-1))) +
    theme_classic() + 
    theme(axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 14), 
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14))
}


plotFigCIs <- function(res, ci, col) {
  setsToPlot <- c("0-1-2-4-8", "0-1-2", "0-2-4", "0-4-8")
  qs <- map(setsToPlot,
            ~ plotCIOrdered(res[[.x]]$fit, ci[[.x]], col) +
              ggtitle(paste(dash2comma(.x), "h")))
  q <- plot_grid(plotlist=qs)
  q
}


getOptTime <- function(k, mu1, mu2, mu3) {
  f <- function(x) {
    pars <- list(d = 1, t = x, mu1 = mu1, mu2 = mu2, mu3 = mu3, k = k)
    I_RTC(pars)
  }
  optimize(f, c(1e-5,20), maximum = TRUE)
}

                       
getIvsTime <- function(pars, times) {
  Is <- map_dbl(times, function(x) {
    pars$t <- x
    I_RTC(pars)
  })
  data.frame(time = times, I = Is)
}


## call


res <- map(dataSets[1:2], function(whichData) {
    rdsDir <- file.path(pulseDir, whichData, fsep=.Platform$file.sep)
    res <- map(dir(rdsDir, paste(modelStr, "fit-0", sep="")), ~readRDS(file.path(rdsDir, .x)))
    timeSets <- extractTimesFromNames(dir(rdsDir, paste(modelStr, "fit-0", sep="")))
    names(res) <- timeSets
    q <- fiForGenesBSA(res[[tSet]])
    ggsave(filename=file.path(figDir, paste(basename, whichData, ".pdf", sep="")),
           plot = q,
           width = 180,
           #height = 210,
           units = "mm",
           dpi = 1200)
    q
})

p1 <- res[[1]]

res <- map(dataSets[3:5], function(whichData) {
    rdsDir <- file.path(pulseDir, whichData, fsep=.Platform$file.sep)
    res <- map(dir(rdsDir, paste(modelStr, "fit-0", sep="")), ~readRDS(file.path(rdsDir, .x)))
    timeSets <- extractTimesFromNames(dir(rdsDir, paste(modelStr, "fit-0", sep="")))
    names(res) <- timeSets
    q <- fiForGenesRTC(res[[tSet]])
    ggsave(filename=file.path(figDir, paste(basename, whichData, ".pdf", sep="")),
           plot = q,
           width = 180,
           #height = 210,
           units = "mm",
           dpi = 1200)
    q
})

p2 <- res[[1]]

# main

g <- plot_grid(p1, p2, labels = "auto", nrow = 2)
ggsave(filename=file.path(figDir, "figure5.pdf", fsep=.Platform$file.sep),
       plot = g,
       width = 180,
       #height = 210,
       units = "mm",
       dpi = 1200)
       

# CIs pulseR

res <- map(dataSets, function(whichData) {
    rdsDir <- file.path(pulseDir, whichData, fsep=.Platform$file.sep)
    
    res <- map(dir(rdsDir, paste(modelStr, "fit-0", sep="")), ~readRDS(file.path(rdsDir, .x)))
    timeSets <- extractTimesFromNames(dir(rdsDir, paste(modelStr, "fit-0", sep="")))
    names(res) <- timeSets
    ci <- map(dir(rdsDir, paste(modelStr, "cis-0", sep="")), ~readRDS(file.path(rdsDir, .x)))
    names(ci) <- timeSets
    
    qRelCIOrdered <- plotFigCIs(res, ci, mpal[[whichData]])
    ggsave(
        file.path(figDir, paste(basename, whichData, ".pdf", sep="")),
        qRelCIOrdered,
        width = 200,
        #height = 160,
        units = "mm",
        dpi = 1200,)
    ggsave(
        file.path(figDir, paste(basename, whichData, ".png", sep="")),
        qRelCIOrdered,
        width = 200,
        #height = 130,
        units = "mm",
        dpi = 1200)
})
    

res <- map(dataSets[grep("slam|tuc|tls", dataSets)], function(whichData) {
    rdsDir <- file.path(pulseDir, whichData, fsep=.Platform$file.sep)
    
    res <- map(dir(rdsDir, paste(modelStr, "fit-0", sep="")), ~readRDS(file.path(rdsDir, .x)))
    timeSets <- extractTimesFromNames(dir(rdsDir, paste(modelStr, "fit-0", sep="")))
    names(res) <- timeSets
    ci <- map(dir(rdsDir, paste(modelStr, "cis-0", sep="")), ~readRDS(file.path(rdsDir, .x)))
    names(ci) <- timeSets
    
    # use tSet
    allPointsFit <- res[["0-1-2-4-8"]]
    allPointsCIs <- ci[["0-1-2-4-8"]]
    
    optRatio <- getOptTime(allPointsFit$fit$size,
                       mu1 = exp(median(allPointsFit$fit$mu1)) + exp(median(allPointsFit$fit$mu2)) + exp(median(allPointsFit$fit$mu3)),
                       mu2 = exp(median(allPointsFit$fit$mu2)),
                       mu3 = exp(median(allPointsFit$fit$mu3)))
                       
    x <- getIvsTime(
        list(d = 1,
            mu1 = exp(median(allPointsFit$fit$mu1)) + exp(median(allPointsFit$fit$mu2)) + exp(median(allPointsFit$fit$mu3)),
            mu2 = exp(median(allPointsFit$fit$mu2)),
            mu3 = exp(median(allPointsFit$fit$mu3)),
            k = allPointsFit$fit$size),
        seq(.01,10,.01))
    
    qInformationVsTime <- qplot(
        data = x,
        x = time,
        y = I,
        geom = "line",
        log='xy',
        xlab = expression(paste(t, "/", tau)),
        ylab = expression(paste(I[paste(delta, delta)] %.% delta^2))) +
        theme(axis.title = element_text(family = "serif"))

    qInformationVsTime <-
        qInformationVsTime +
        annotate("point",
            x = optRatio$maximum,
            y = optRatio$objective, shape = 1, size = 3,
        colour = '#3182BD') +
        annotate("segment",
            x = optRatio$maximum,
            xend = optRatio$maximum,
            y = 0,
            yend = optRatio$objective,  size = .5,
            colour = '#3182BD') +
        scale_x_continuous(
            trans = "log",
            breaks = c(.01,.05,optRatio$maximum, 1, 5, 10),
            labels = c("0.01","0.05", format(optRatio$maximum, digits=2), "1", "5", "10")) +
        ggtitle(toupper(whichData)) +
        theme_classic() + 
                    theme(axis.text.y = element_text(size = 16), 
                        axis.title.y = element_text(size = 20),
                        axis.text.x = element_text(size = 16), 
                        axis.title.x = element_text(size = 20))
        qInformationVsTime
})


g <- plot_grid(res[[1]], 
               res[[2]] + 
                    theme(axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title.y = element_blank() ), 
               res[[3]] + 
                    theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank() ),
                   nrow = 1, ncol = 3,
                   labels = "auto",
                   align = "v")
                        
ggsave(filename=file.path(figDir, "optt.pdf", fsep=.Platform$file.sep),
       plot = g,
       width = 350,
       #height = 210,
       units = "mm",
       dpi = 1200)
       
