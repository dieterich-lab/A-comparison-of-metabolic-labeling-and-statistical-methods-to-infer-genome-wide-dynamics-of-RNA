
## various 

extractTimesFromNames <- function(files) {
  timeSets <- stringr::str_match_all(files, "-([0-9\\.].+)\\.rds") %>%
    map(2)
  timeSets
}


naToVal <- function(x, val) {
  x[is.na(x)] <- val
  x
}


naToBorders <- function(x, vals) {
  for(i in 1:2)
    x[,i] <- naToVal(x[,i], vals[i])
  x
}


widths <- function(x)
    x[,2] - x[,1]
    
    
relWidth <- function(x,d)
    widths(x)/d

    
points2label <- function(times)
  paste(unique(times), collapse = "-")
  
  
dash2comma <- function(x)
  str_replace_all(x, "-", ",")
  

lappend <- function (lst, ...){
    lst <- c(lst, list(...))
    return(lst)
  }


# FIM terms for BSA (standard model) and SLAM/TLS/TUC model
# labelled dd FIM term
I_labelledBSA <- function(pars, NB = TRUE) {
    eval(expression({
        I <- t^2 * exp(-2*d*t) * mu/(1 - exp(-d * t))
        if (NB) { 
            I <- I * 1/(1 + mu/k * (1 - exp(-d * t))) 
        }
        I
    }),
    envir = pars)
}


# unlabelled dd FIM term
I_unlabelledBSA <- function(pars, NB = TRUE) {
    eval(expression({
        I <- t^2 * mu * exp(-d*t)
        if (NB) { 
            I <- I /(1 + mu/k *  exp(-d*t)) 
        }
        I
    }), 
    pars)
}


# labelled dd FIM term - for RTC model
I_labelledRTC <- function(pars, NB = TRUE) {
    eval(expression({
        I <- t^2 * exp(-2 * d * t) * mu3^2/(mu2 + mu3 * (1 - exp(-d * t)))
        if (NB) { 
            I <- I * 1/(1 + mu2/k  + mu3/k * (1 - exp(-d * t))) 
        }
        I
    }),
    envir = pars)
}


# unlabelled dd FIM term - for RTC model
I_unlabelledRTC <- function(pars, NB = TRUE) {
    eval(expression({
        I <- t^2 * exp(-2 * d * t) * mu3^2/(mu1 + mu3 * exp(-d * t))
        if (NB) { 
            I <- I * 1/(1 + mu1/k + mu3/k * exp(-d * t)) 
        }
        I
    }), 
    pars)
}


I_RTC <- function(pars) {
  I_labelledRTC(pars, NB = TRUE) + I_unlabelledRTC(pars, NB = TRUE)
}


## plotting

shadesOfGrey <- colorRampPalette(c("white", "grey100", "grey90", "grey80", "grey70", "grey60", 
                                   "grey50", "grey40", "grey30", "grey20", "grey10", "grey0"))
                                   
                                   
myHeatmap <- function(cormat, labels, llim=0.6) {
    # melt fo plotting...
	melted_cormat <- reshape2::melt(cormat, na.rm = TRUE)
	# manually set column order 
    col_order <- rev(levels(melted_cormat$Var2))
    melted_cormat$Var1 <- factor(melted_cormat$Var1, levels = col_order)
    melted_cormat$Var2 <- factor(melted_cormat$Var2, levels = col_order)
	# and create a ggheatmap
	ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
        geom_tile(color = "white") +
        scale_fill_gradient(low = "#ece7f2", high = "#3182BD", limit = c(llim,1), space = "Lab", 
            name="Pearson\nCorrelation",na.value='white',guide=F) +
        theme_minimal() + 
        theme(text=element_text(size=12),
              axis.text.y = element_text(size=12),
              axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12, hjust = 1)) +
		scale_x_discrete(labels=labels) + scale_y_discrete(labels=labels) +
        coord_fixed()
	# add text (corr values)
	ggheatmap + 
        geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
	          axis.ticks = element_blank())
}                                   
                                   
                                   
myLine <- function(x, y, ...){
    par(new = TRUE)
    smoothScatter(x, y, 
                  colramp = shadesOfGrey,
                  xaxt='n', yaxt = 'n',
                  ...)
    abline(a = 0,b = 1, ...)
}


# Viewport function
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}


log10.axis <- function(side, at, ...) {
    at.minor <- log10(outer(1:9, 10^(min(at):max(at))))
    axis(side=side, at=at.minor, labels=NA, tcl=par("tcl")*0.5, ...)
}


lower.panel <- function(x, y, i, j, ...){
    # if comparing to ERCC, use residual SE from linear fit
    # ERCC always need to be on first column!
    # else use RMSD
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    if (j==1) {
        f <- lm(y ~ x)
        rse <- round(sqrt(deviance(f)/f$df.residual), 2)
        rse <- format(rse, digits=2, nsmall=2)
        txt <- paste0("RSE\n", rse)
    } else {
        rmsd <- round(sqrt(mean((y - x)^2)), 2)
        rmsd <- format(rmsd, digits=2, nsmall=2)
        txt <- paste0("RMSD\n", rmsd)
    }
    cex.cor <- .8/strwidth(txt)
    text(grconvertX(0.5,"npc"), grconvertY(0.5, "npc"), txt, cex = cex.cor)
}


# passing indices to lower panels
# adapted from R function pairs
myPairs <-

function (x, labels, panel = points, ...,

          horInd = 1:nc, verInd = 1:nc,

          lower.panel = panel, upper.panel = panel,

          diag.panel = NULL, text.panel = textPanel,

          label.pos = 0.5 + has.diag/3, line.main = 3,

          cex.labels = NULL, font.labels = 1,

          row1attop = TRUE, gap = 1, log = "",

          horOdd = !row1attop, verOdd = !row1attop)

{

    if(doText <- missing(text.panel) || is.function(text.panel))

	textPanel <-

	    function(x = 0.5, y = 0.5, txt, cex, font)

		text(x, y, txt, cex = cex, font = font)

 

    localAxis <- function(side, x, y, xpd, bg, col=NULL, main, oma, ...) {

      ## Explicitly ignore any color argument passed in as

      ## it was most likely meant for the data points and

      ## not for the axis.

        xpd <- NA

        if(side %% 2L == 1L && xl[j]) xpd <- FALSE

        if(side %% 2L == 0L && yl[i]) xpd <- FALSE

        if(side %% 2L == 1L) Axis(x, side = side, xpd = xpd, ...)

        else Axis(y, side = side, xpd = xpd, ...)

    }

 

    localPlot <- function(..., main, oma, font.main, cex.main) plot(...)

    localLowerPanel <- function(..., main, oma, font.main, cex.main)

        lower.panel(...)

    localUpperPanel <- function(..., main, oma, font.main, cex.main)

        upper.panel(...)

 

    localDiagPanel <- function(..., main, oma, font.main, cex.main)

        diag.panel(...)

 

    dots <- list(...); nmdots <- names(dots)

    if (!is.matrix(x)) {

        x <- as.data.frame(x)

        for(i in seq_along(names(x))) {

            if(is.factor(x[[i]]) || is.logical(x[[i]]))

               x[[i]] <- as.numeric(x[[i]])

            if(!is.numeric(unclass(x[[i]])))

                stop("non-numeric argument to 'pairs'")

        }

    } else if (!is.numeric(x)) stop("non-numeric argument to 'pairs'")

    panel <- match.fun(panel)

    if((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))

        lower.panel <- match.fun(lower.panel)

    if((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))

        upper.panel <- match.fun(upper.panel)

    if((has.diag  <- !is.null( diag.panel)) && !missing( diag.panel))

        diag.panel <- match.fun( diag.panel)

 

    if(row1attop) {

        tmp <- lower.panel; lower.panel <- upper.panel; upper.panel <- tmp

        tmp <- has.lower; has.lower <- has.upper; has.upper <- tmp

    }

 

    nc <- ncol(x)

    if (nc < 2L) stop("only one column in the argument to 'pairs'")

    if(!all(1L <= horInd & horInd <= nc))

        stop("invalid argument 'horInd'")

    if(!all(1L <= verInd & verInd <= nc))

        stop("invalid argument 'verInd'")

    if(doText) {

	if (missing(labels)) {

	    labels <- colnames(x)

	    if (is.null(labels)) labels <- paste("var", 1L:nc)

	}

	else if(is.null(labels)) doText <- FALSE

    }

    oma  <- if("oma"  %in% nmdots) dots$oma

    main <- if("main" %in% nmdots) dots$main

    if (is.null(oma))

	oma <- c(4, 4, if(!is.null(main)) 6 else 4, 4)

    opar <- par(mfcol = c(length(horInd), length(verInd)),

                mar = rep.int(gap/2, 4), oma = oma)

    on.exit(par(opar))

    dev.hold(); on.exit(dev.flush(), add = TRUE)

 

    xl <- yl <- logical(nc)

    if (is.numeric(log)) xl[log] <- yl[log] <- TRUE

    else {xl[] <- grepl("x", log); yl[] <- grepl("y", log)}

    ni <- length(iSet <- if(row1attop) horInd else rev(horInd))

    nj <- length(jSet <- verInd)

    for(j in jSet)

        for(i in iSet) {

            l <- paste0(if(xl[j]) "x" else "",

                        if(yl[i]) "y" else "")

            localPlot(x[, j], x[, i], xlab = "", ylab = "",

                      axes = FALSE, type = "n", ..., log = l)

            if(i == j || (i < j && has.lower) || (i > j && has.upper) ) {

                box()

                j.odd <- (match(j, jSet) + horOdd) %% 2L

                i.odd <- (match(i, iSet) + verOdd) %% 2L

                if(i == iSet[1L] && (!j.odd || !has.upper || !has.lower))

                    localAxis(3L, x[, j], x[, i], ...)

                if(i == iSet[ni] && ( j.odd || !has.upper || !has.lower))

                    localAxis(1L, x[, j], x[, i], ...)

                if(j == jSet[1L] && (!i.odd || !has.upper || !has.lower))

                    localAxis(2L, x[, j], x[, i], ...)

                if(j == jSet[nj] && ( i.odd || !has.upper || !has.lower))

                    localAxis(4L, x[, j], x[, i], ...)

                mfg <- par("mfg")

                if(i == j) {

                    if (has.diag) localDiagPanel(as.vector(x[, i]), ...)

		    if (doText) {

                        par(usr = c(0, 1, 0, 1))

                        if(is.null(cex.labels)) {

                            l.wid <- strwidth(labels, "user")

                            cex.labels <- max(0.8, min(2, .9 / max(l.wid)))

                        }

                        xlp <- if(xl[i]) 10^0.5 else 0.5

                        ylp <- if(yl[j]) 10^label.pos else label.pos

                        text.panel(xlp, ylp, labels[i],

                                   cex = cex.labels, font = font.labels)

                    }

                } else if(i < j)

                    localLowerPanel(as.vector(x[, j]), as.vector(x[, i]), ...)

                else

                    localUpperPanel(as.vector(x[, j]), as.vector(x[, i]), i, j, ...)

                if (any(par("mfg") != mfg))

                    stop("the 'panel' function made a new plot")

            }

            else par(new = FALSE)

        }

    if (!is.null(main)) {

        font.main <- if("font.main" %in% nmdots) dots$font.main else par("font.main")

        cex.main  <- if("cex.main"  %in% nmdots) dots$cex.main  else par("cex.main")

        mtext(main, 3, line.main, outer=TRUE, at = 0.5, cex = cex.main, font = font.main)

    }

    invisible(NULL)

}


# customized histogram from the Cellplot package
# https://github.com/dieterich-lab/CellPlot
myHistogram = function(framelist, go.alpha=0.05, go.alpha.term="pvalCutOff", 
                        min.sig=1, min.genes=10, max.genes=100, bar.scale=NULL,
                        reorder=T, main="GO enrichment", show.go.id=FALSE, prefix="", 
                        lab.cex=1, axis.cex=1, group.cex=NULL, go.selection=NULL, 
                        term.selection=NULL, main.cex=1, mycolor1=mycolor1, mycolor2=mycolor2, fixed.lim=NULL) {
  
    grep.multi = function(xvec, target) {
      acc = vector(mode="numeric")
      for (i in 1:length(xvec)) {
        acc = c(acc, grep(xvec[i], target))
      }
      return(unique(acc))
    }
  
  # same pattern, we don't care about the names, we just use 'LogEnriched', and put it in 'Upregulated'
  for (i in 1:length(framelist)) {
    framelist[[i]]$Upregulated = framelist[[i]]$Downregulated = 0 #framelist[[i]]$SignificantValues = 0
    framelist[[i]]$p.up = framelist[[i]]$p.down = NA
    for (j in 1:nrow(framelist[[i]])) {
      framelist[[i]]$Upregulated[j] = framelist[[i]]$LogEnriched[j]
    }
  }
  
  # preprocess into one combined data frame
    if (!is.null(go.selection)) {
      for (f in 1:length(framelist)) {
        framelist[[f]] = framelist[[f]][which(framelist[[f]]$GO.ID %in% go.selection),]
      }
    }
    if (!is.null(term.selection)) {
      for (f in 1:length(framelist)) {
        framelist[[f]] = framelist[[f]][which(framelist[[f]]$Term %in% term.selection),]
      }
    }
    if(is.null(names(framelist))) { names(framelist) = paste0("group_",1:length(framelist) ) }
    framelist = lapply(framelist, function(x) { x$idterm = paste(x$GO.ID, x$Term); return(x) })
  
  # Integrate the topGO frames into one combined data.frame
    allrows = vector()
    for (i in 1:length(framelist)) { 
      sel = which(!is.na(framelist[[i]][,"idterm"]))
      framelist[[i]] = framelist[[i]][sel,]
      rownames(framelist[[i]]) = framelist[[i]][,"idterm"]
      allrows = c(allrows, rownames(framelist[[i]]))
    }
    allrows = unique(allrows)
    
    intframe = data.frame(row.names=allrows)
    outcolnames = vector()
    for (j in 1:length(framelist)) {
      tmp = data.frame( framelist[[j]] )
      colnames(tmp) = colnames( framelist[[j]] )
      rownames(tmp) = rownames(framelist[[j]])
      outcolnames = c( outcolnames, paste( names(framelist[j]), colnames(tmp), sep="." ) )
      intframe = data.frame( cbind(intframe, tmp[rownames(intframe),]) )
      rownames(intframe) = allrows
    }
    colnames(intframe) = outcolnames
  
  # select terms if specified
    if (!is.null(term.selection)) {
      termorder = gsub( "GO:[0-9]+ ","", rownames(intframe) )
      termorder = match(term.selection, termorder)
      intframe = intframe[rev(termorder),]
    }
  
  # apply filters
    filterframe = (intframe[,grep(paste0("\\.", go.alpha.term),colnames(intframe))] <= go.alpha) *  ( intframe[,grep("\\.Significant$",colnames(intframe))] >= min.genes )
    intframe = intframe[ apply(filterframe,1,function(x) sum(x, na.rm = T)) >= min.sig, ]
    intframe = intframe[ which(apply( intframe[,grep("\\.Significant$",colnames(intframe))], 1, function(x) all( x <= max.genes, na.rm = T) ) ), ]
    # select terms with at least min.sig samples
    
  # cluster by significance status (p[enrich]<=alpha) across groups (other options?)
    if (reorder) { 
      sf = intframe[,grep(paste0("\\.", go.alpha.term),colnames(intframe))]
      sf = sf <= go.alpha
      sf[is.na(sf)] = FALSE
      sfo = hclust(dist(sf))$order
      intframe = intframe[sfo,] 
    }
  
  # limit for axis
    max.lim <- max(intframe[,grep.multi(paste0(prefix, c("LogEnriched")),colnames(intframe))], na.rm=T)
    print(max.lim)
    if (!is.null(fixed.lim)){
        lim <- fixed.lim
    } else {
        lim <- max.lim
    }
  
  # Layout
    dv = dev.size(units = "in")
    m = m.old = par("mai")
    m[c(1,3)] = sapply( m[c(1,3)]/dv[2], function(x) max(x, 0.01) )
    m[c(2,4)] = sapply( m[c(2,4)]/dv[1], function(x) max(x, 0.01) )
    nc = length(framelist) * 2      # n.o. columns, just a shortcut variable
    thi = 0.6; th = thi / dv[2]                # desired label area height in inches (then scaled)
    bhi = 0.18; bh = (bhi * ifelse(is.null(bar.scale),1,bar.scale))  / dv[2]       # height of one bar if fixed in inches
    ta = 2 / dv[1]    # label area in inches (then scaled)
    bot = 0           # bottom margin of bar plots
    top = 0.4         # top margin of bar plots (need to accomodate the axes)
    gap = 0.05         # space to the left/right of bars
  
    if(!is.null(bar.scale)) { 
      neededspace = bhi * nrow(intframe) + sum( m.old[c(1,3)] ) + thi
      if ( neededspace > dv[2] ) { stop(paste0("figure region too small (set device height to above ",round(neededspace,digits = 2)," in)")) }
    }
  
    layoutmatrix = rbind( c( nc/2+nc+5, rep( nc/2+nc+5, nc+1) ),                # m[3] # top margin  
                          c( nc/2+nc+4, sort( rep( 1:(nc/2), 2) ), nc/2+1 ),    # th
                          c( nc/2+nc+4, (nc/2+2):(nc/2+nc+1), nc/2+nc+2 ),      # bararea
                          c( nc/2+nc+3, rep( nc/2+nc+3, nc+1) ) )               # m[1] # bottom margin
  
    barareaheight = ifelse( is.null(bar.scale), 
                            1-m[3]-m[1]-th, 
                            min( (top+bot)/dv[2]+bh*nrow(intframe), 1-m[3]-m[1]-th ) )
    barplotwidth = (1-m[2]-m[4]-ta)/nc
  
    if ( any( c(barareaheight, barplotwidth-gap/dv[1])  <= 0 ) ) { 
        par(mai=m.old, xpd = F); layout(1); stop("figure region too small") 
    }
  
  layout(layoutmatrix, 
           heights = c( m[3],th,barareaheight,1-m[3]-th-barareaheight),
           widths = c( m[2], rep( barplotwidth, nc ), m[4]+ta ) )

  # group labels and title
    par(mai=c(0,0,0,0), xpd = T)
    for (i in 1:length(framelist) ) {
      plot.new()
      text(0.5,0.5,unlist(strsplit(names(framelist)[i], " "))[1], cex = ifelse(is.null(group.cex), 2.5, group.cex), srt=0)
    }
    plot.new()
    text(0.5,0.5,main, cex=main.cex)
  
  # Cycling through the groups and plotting each in turn
  # some lines of codes are left untouched, but ignored in the end for the plot
  # we just leave them for now...
                       
  for (i in 1:length(framelist) ) {
    #tmp = intframe[,grep(names(framelist[i]), colnames(intframe))]
    # HARD CODED for our purpose, else this is really unstable piece of code!
    tmp = intframe[,grep(paste(names(framelist[i]), '.[A-Z,a-z]', sep=''), colnames(intframe))]
    ts = tmp[,grep(paste0("\\.", go.alpha.term),colnames(tmp))]
    ucol = rep("darkgrey",length(ts))
    dcol = rep("grey",length(ts))
    mycolor <- mycolor1
    if (i > length(framelist)/2) { mycolor <- mycolor2 }
    ucol[which(ts <= go.alpha)] = mycolor
    dcol[which(ts <= go.alpha)] = "deepskyblue2"
    
    # left-pointed bars / downregulated
    # nothing here, just keep the lines ...
      par(mai=c(bot,0,top,0))
      
      vline = barplot(-tmp[,grep(paste0(prefix,"Downregulated"),colnames(tmp))], horiz = T, 
                      xlim = c(0,0), border=NA, col = dcol, axes=F)       # getting the bar positions
      abline(h=vline, col="grey", lty=2)

    # right-pointed bars / upregulated, and middle lines
      par(mai=c(bot,0,top,gap))
      barplot(tmp[,grep(paste0(prefix,"Upregulated"),colnames(tmp))], horiz = T, 
              xlim = c(0,lim), border=NA, col=ucol, axes=F)
      abline(h=vline, col="grey", lty=2)
      barplot(tmp[,grep(paste0(prefix,"Upregulated"),colnames(tmp))], horiz = T, 
              xlim = c(0,lim), border=NA, col=ucol, axes=F, add=T)
      axis(3, at=floor( seq(0,lim,length.out = 2)), cex.axis=axis.cex, las=3 )
      lines(c(0,0),c( par('usr')[4], -0.5*vline[1] ) )
    
  }
  
  
  # Row labels
    par(mai=c(bot,0,top,0))
    barplot( rep(NA,length(vline)), axes = F, horiz=T, xlim = c(0.1,1) ) # dummy plot to get the scale -- any better way?
    text(0.1, vline, pos = 4, cex=lab.cex,
         labels = if (show.go.id) { rownames(intframe) } else { gsub("GO:[0-9]+ ","",rownames(intframe)) })

  # bottom margin and left margin placeholder plots
    par(mar=c(0,0,0,0) )
    plot.new()
    plot.new()
    plot.new()
  
  # Reset graphical params
    par(mai=m.old, xpd = F)
    layout(1)
  
}


# custom
statBoxData <- function(y, upper_limit = -0.15) { # if top -> upper_limit = max(df$d) * 1.15
  return( 
    data.frame(
      y = upper_limit,
      label = length(y)
    )
  )
}


# custom 
myBoxPlotsPlot <- function(bdat, att, att_color){
  yaxis <- as.character(att)
  if (grepl("^d.", yaxis)) {
      ylab <- expression(paste("Median ", delta, ", ", h^-1))
  } else {
      ylab <- "Median CI width"
  }   
  col <- match(att, colnames(bdat))
  colnames(bdat)[col] <- "attribute"
  upper_xlim <- as.numeric((max(bdat$x) + 1))
  plot_lims <- as.numeric(0:upper_xlim)
  bdat$x <- as.factor(bdat$x)
  boxplots <- ggplotGrob(ggplot()
                         + theme_bw() + ylab(ylab)
                         + scale_x_discrete(limits = plot_lims, expand = c(0,0))
                         + theme(plot.margin = unit(c(-0.7,0,0,0), "cm"),
                                 axis.title.y = element_text(vjust = 0.2, size = 22, angle = 0),
                                 axis.ticks.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 panel.border = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 panel.grid.major = element_blank(),
                                 axis.title.x = element_blank(),
                                 axis.text.y = element_text(size = 16))
                         + geom_boxplot(data = bdat, aes_string(x="x", y="attribute"),
                                        fill = att_color, colour = "gray80"))
    
  
  return(boxplots)
}


# Generates UpSet plot with boxplots representing distributions of attributes
myBaseBoxPlot <- function(box_plot, position, size_plot_height, Main_bar_plot, Matrix_plot,
                        Size_plot, hratios, set_metadata, set_metadata_plots, newpage){

  if(length(box_plot) > 2){
    return(warning("UpSet can only show 2 box plots at a time"))
  }
  if(is.null(position) || position == tolower("bottom")){
    bar_top <- 1
    matrix_bottom <- 100
    att_top <- 101
    att_bottom <- 130
    if(length(box_plot) == 2){
      att_top <- 105
      att_bottom <- 120
      gridrow <- 145
    }
  }
  if((!is.null(position)) && (position != tolower("bottom"))){
    if(length(box_plot)==1){
      size_plot_height <- (size_plot_height + 35)
      bar_top <- 36
      matrix_bottom <- 135
      att_top <- 10
      att_bottom <- 35
    }
    else if(length(box_plot) == 2){
      size_plot_height <- (size_plot_height + 50)
      bar_top <- 51
      matrix_bottom <- 150
      att_top <- 15
      att_bottom <- 30
      gridrow <- 150
    }
  }
  if(is.null(set_metadata)){
    matrix_and_mainbar_right <- 100
    matrix_and_mainbar_left <- 21
    size_bar_right <- 20
    size_bar_left <- 1
  }
  else if(!is.null(set_metadata)){
    matrix_and_mainbar_right <- set_metadata$ncols + 100
    matrix_and_mainbar_left <- set_metadata$ncols + 21
    size_bar_right <- set_metadata$ncols + 20
    size_bar_left <- set_metadata$ncols + 1
    metadata_right <- set_metadata$ncols
    metadata_left <- 1
  }
  if (newpage) {
    grid.newpage()
  }
  if(length(box_plot) == 1){
    pushViewport(viewport(layout = grid.layout(135,matrix_and_mainbar_right)))
  }
  else if(length(box_plot) == 2){
    pushViewport(viewport(layout = grid.layout(gridrow,matrix_and_mainbar_right)))
  }
  vp = vplayout(bar_top:matrix_bottom, matrix_and_mainbar_left:matrix_and_mainbar_right)
  pushViewport(vp)
  grid.draw(arrangeGrob(Main_bar_plot, Matrix_plot, heights = hratios))
  popViewport()
  vp = vplayout(size_plot_height:matrix_bottom, size_bar_left:size_bar_right)
  pushViewport(vp)
  grid.draw(arrangeGrob(Size_plot))
  popViewport()
  if(!is.null(set_metadata)){
    for(i in 1:length(set_metadata_plots)){
      if(i != 1){
        metadata_left <- 1+metadata_right
        metadata_right <- metadata_right + set_metadata$plots[[i]]$assign
      }
      else{
        metadata_left <- 1
        metadata_right <- set_metadata$plots[[i]]$assign
      }

      vp = vplayout(size_plot_height:matrix_bottom, metadata_left:metadata_right)
      pushViewport(vp)
      grid.draw(arrangeGrob(set_metadata_plots[[i]]))
      popViewport()
    }
  }
  vp = vplayout((att_top - 5):att_bottom, (matrix_and_mainbar_left+2):(matrix_and_mainbar_right-2))
  pushViewport(vp)
  grid.draw(arrangeGrob(box_plot[[1]]))
  popViewport()
  if(length(box_plot) == 2){
    vp = vplayout((att_bottom + 5):(att_bottom + 25), (matrix_and_mainbar_left+2):(matrix_and_mainbar_right-2))
    pushViewport(vp)
    grid.draw(arrangeGrob(box_plot[[2]]))
    popViewport()
  }
}


# adapted from https://github.com/dieterich-lab/DesignMetabolicRNAlabeling
# plot scatter of counts vs. predicted (pulseR) for standard model
plotBySample <- function(fit, pd, conditions) {
  library(tidyverse)
  library(pulseR)
  pr <- predictExpression(fit, pd)
  # adjust names
  conditions[grep("flow", conditions$fraction),]$fraction <- "Supernatant"
  conditions[grep("pull", conditions$fraction),]$fraction <- "Eluate"
  conditions[grep("total", conditions$fraction),]$fraction <- "Input"
  counts <- as.data.frame(pd$counts)
  counts <- counts %>% rownames_to_column("gene_id") %>%
    gather(sample, count, -gene_id) %>%
    inner_join(conditions)
  predictions <- data.frame(pr$predictions) %>%
    `row.names<-`(rownames(pd$counts)) %>%
    set_names(colnames(pd$counts)) %>%
    rownames_to_column("gene_id") %>%
    gather(sample, count, -gene_id) %>%
    inner_join(conditions)
  labs <- unite(conditions, "S", fraction, rep, time, sep = " ")$S
  names(labs) <- conditions$sample
  counts$type <- "Measured"
  predictions$type <- "Estimated"
  counts$mu <- predictions$mu <- exp(fit$mu)
  counts$d <- predictions$d <- fit$d
  x <- rbind(counts, predictions)
  x$newCol <- paste(x$fraction, x$rep, x$time, sep="")
  x %>% spread(type, count) %>%
    ggplot(data = .) +
    geom_point(aes(
      x = Estimated,
      y = Measured,
      colour = d
    ), size = .2) +
    geom_abline(slope = 1,
                intercept = 0,
                colour = "black") +
    facet_wrap( ~ newCol) + 
    scale_y_log10() +
    scale_x_log10() +
    scico::scale_colour_scico(palette = 'grayC', direction=-1) +
    theme(strip.background = element_rect(fill = NA),
          strip.text.x = element_text(margin = margin(.2,0,.2,0, "cm"))) +
    guides(x = guide_axis(angle = 45)) +
    theme_minimal() + labs(color=substitute(paste(delta, " ", h^-1)))
  
}

