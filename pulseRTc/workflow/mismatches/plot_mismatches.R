#!/biosw/R/3.5.1/bin/Rscript

# plot mismatches 
# we follow same nomenclature/structure as GRAND-SLAM for comparison

# Usage: ./plot_mismatches.R <slam/tls/tuc> <1/2>
# 1: <slam/tls/tuc>
# 2: <1/2> 1: used mismatches (default), 2: all mismatches (does not generate mismatchpos.pdf) 


library(ggplot2)
library(plyr)
library(purrr)
library(data.table)

# get params and options
args <- commandArgs(trailingOnly=TRUE)
which <- TRUE
if (length(args)<1) {
  stop("USAGE: plot_mismatches.R <slam/tls/tuc>\n", call.=FALSE)
} else {
  if (length(args)>1 & as.integer(args[2])==2) {
    which <- FALSE
  }
}

loc <- here::here()
prefix <- 'pulseRTc/workflow/mismatches'
prefix <- file.path(loc, prefix, args[1], fsep=.Platform$file.sep)
print(paste("Processing ", prefix, " ...", sep=""))

outdir <- prefix
print(paste("Writing to ", outdir, " ...", sep=""))

basename <- paste(toupper(args[1]), "mismatches.pdf", sep='.')

# first construct full dataframe
mismatches <- map(list.files(prefix, pattern = 'mismatches.tab.gz$', full.names = T), read.table, head = TRUE)
names <- vapply(strsplit(list.files(prefix, pattern = 'mismatches.tab.gz$'), ".", fixed = TRUE), "[", "", 1)
names(mismatches) <- names
for(i in seq_along(mismatches)) {
    mismatches[[i]]$Condition <- names[i]
}
all.mismatches <- do.call("rbind", mismatches)

if (which) {
    used <- map(list.files(prefix, pattern = 'mismatches-used.tab.gz$', full.names = T), read.table, head = TRUE)
    names <- vapply(strsplit(list.files(prefix, pattern = 'mismatches-used.tab.gz$'), ".", fixed = TRUE), "[", "", 1)
    names(used) <- names
    for(i in seq_along(used)) {
        used[[i]]$Condition <- names[i]
        used[[i]] <- used[[i]] %>% dplyr::rename(Used=Mismatches)
    }
    all.used <- do.call("rbind", used)


    fields = colnames(all.used)
    fields = fields[fields!='Used']

    all.mismatches <- merge(all.mismatches, all.used, by=fields, all = T)
    # reorder
    o <- order(all.mismatches$Condition, all.mismatches$Orientation, all.mismatches$Genomic, all.mismatches$Read)
    all.mismatches <- all.mismatches[o,]
    
    all.mismatches <- all.mismatches %>%
                        dplyr::mutate(Used = dplyr::coalesce(Used ,Mismatches)) %>%
                        dplyr::select(-Mismatches) %>% dplyr::rename(Mismatches=Used)
    basename <- paste(toupper(args[1]), "mismatches-used.pdf", sep='.')
}

# prep dataframe
all.mismatches$Mismatch <- paste0(all.mismatches$Genomic,"->",all.mismatches$Read)
all.mismatches$Mismatch <- factor(all.mismatches$Mismatch,levels=unique(all.mismatches$Mismatch))
all.mismatches$Rate <- all.mismatches$Mismatches/all.mismatches$Coverage
all.mismatches$se <- sqrt(all.mismatches$Rate*(1-all.mismatches$Rate)/all.mismatches$Coverage)
all.mismatches$Condition <- factor(all.mismatches$Condition,as.character(unique(all.mismatches$Condition)))

pdf(file.path(outdir, basename, fsep=.Platform$file.sep), width=4+length(unique(all.mismatches$Condition)), height=7)
# mismatch statistics for all first and second reads (whether sense or antisense reads)
# currently Any (exonic, intronic)
for (category in unique(all.mismatches$Category)) {
	for (or in unique(all.mismatches$Orientation)) {
		if (sum(all.mismatches$Category==category & all.mismatches$Orientation==or)>0) {
            print(ggplot(all.mismatches[all.mismatches$Category==category & all.mismatches$Orientation==or,],
                  aes(Mismatch, Rate, color=Condition)) + 
                  geom_point(position=position_dodge(width=0.7)) + 
                  geom_errorbar(aes(ymin=Rate-se, ymax=Rate+se), alpha=0.3, position=position_dodge(width=0.7), width=0) + 
                  theme_bw() + 
                  theme(text=element_text(size=24), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
                  ggtitle(paste(or)))
        }
    }
}
# same, but by mismatch type
for (mm in unique(all.mismatches$Mismatch)) {
	for (or in unique(all.mismatches$Orientation)) {
		if (sum(all.mismatches$Mismatch==mm & all.mismatches$Orientation==or)>0) {
            print(ggplot(all.mismatches[all.mismatches$Mismatch==mm & all.mismatches$Orientation==or,],
            aes(Category, Rate, color=Condition)) + 
            geom_point(position=position_dodge(width=0.7)) + 
            geom_errorbar(aes(ymin=Rate-se,ymax=Rate+se), alpha=0.3, position=position_dodge(width=0.7), width=0) + 
            theme_bw() + 
            theme(text=element_text(size=24), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
            ggtitle(paste(mm, or)))
		}
	}
}	
dev.off()

if (which) {
    # same, we construct full dataframe, but we plot each sample separately
    mismatches <- map(list.files(prefix, pattern = 'mismatchDetails.tab.gz$', full.names = T), read.table, head = TRUE)
    names <- vapply(strsplit(list.files(prefix, pattern = 'mismatchDetails.tab.gz$'), ".", fixed = TRUE), "[", "", 1)
    names(mismatches) <- names
    for(i in seq_along(mismatches)) {
        mismatches[[i]]$Condition <- names[i]
    }
    all.mismatches <- do.call("rbind", mismatches)
    # prep dataframe
    all.mismatches$Mismatch <- paste0(all.mismatches$Genomic,"->",all.mismatches$Read)
    all.mismatches$Mismatch <- factor(all.mismatches$Mismatch,levels=unique(all.mismatches$Mismatch))
    all.mismatches$Rate <- all.mismatches$Mismatches/all.mismatches$Coverage
    all.mismatches$Condition <- factor(all.mismatches$Condition,as.character(unique(all.mismatches$Condition)))
    
    basename <- paste(toupper(args[1]), "mismatchpos.pdf", sep='.')
    pdf(file.path(outdir, basename, fsep=.Platform$file.sep), width=10, height=6)
    for (cond in unique(all.mismatches$Condition)) {
        for (category in unique(all.mismatches$Category)) {
            print(ggplot(all.mismatches[all.mismatches$Category==category & all.mismatches$Condition==cond,],
                aes(Position, Rate, color=Read)) + 
                geom_line() + 
                facet_grid(~Genomic) + 
                theme_bw() + 
                theme(text=element_text(size=24), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
                ggtitle(cond))
        }
    }
    dev.off()
}

